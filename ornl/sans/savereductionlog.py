from datetime import datetime
import h5py
import json
from mantid import __version__ as mantid_version
from mantid.simpleapi import mtd, SaveNexusProcessed
import numpy as np
from ornl import __version__ as sangria_version
from os import environ


def _createnxgroup(parent, name, klass):
    child = parent.create_group(name)
    child.attrs['NX_class'] = klass
    return child


def _savenxnote(nxentry, name, mimetype, file_name, data):
    '''
    http://download.nexusformat.org/doc/html/classes/base_classes/NXnote.html

    :param nxentry:
    :param name:
    :param mimetype:
    :param file_name:
    :param data:
    :return: The generated note
    '''
    nxnote = _createnxgroup(nxentry, name, 'NXnote')

    nxnote.create_dataset(name='file_name', data=[np.string_(file_name)])
    nxnote.create_dataset(name='type', data=[np.string_(mimetype)])
    nxnote.create_dataset(name='data', data=[np.string_(data)])

    return nxnote


def _savepythonscript(nxentry, **kwargs):
    '''
    :param entry: The h5py entry that the note should be in
    :param kwargs: A dict containing ``python`` and ``pythonfile`` arguments.
    See :ref:`savereductionlog` for what they mean.
    '''
    # read the contents of the script
    scriptfile = kwargs.get('pythonfile', '')
    scriptcontents = kwargs.get('python', '')
    if not scriptcontents:
        with open(scriptfile, 'r') as script:
            scriptcontents = script.read()

    return _savenxnote(nxentry, 'reduction_script', 'text/x-python',
                       scriptfile, scriptcontents)


def _savereductionparams(nxentry, **kwargs):
    parameters = kwargs.get('reductionparams', '')
    if not isinstance(parameters, str):
        parameters = json.dumps(parameters)
    return _savenxnote(nxentry, 'reduction_parameters',
                       'application/json', file_name='',
                       data=parameters)


def _savenxprocess(nxentry, program, version):
    nxprocess = _createnxgroup(nxentry, program, 'NXprocess')
    nxprocess.create_dataset(name='program', data=[np.string_(program)])
    nxprocess.create_dataset(name='version', data=[np.string_(version)])


def _savenxlog(nxcollection, property):
    nxlog = _createnxgroup(nxcollection, property.name, 'NXlog')

    try:
        if isinstance(property.value, str):
            value = nxlog.create_dataset(name='value',
                                         data=[np.string_(property.value)])
        elif len(property.value) > 1:
            value = nxlog.create_dataset(name='value', data=property.value)
        else:
            raise RuntimeError('Should never have gotten here')
    except TypeError:
        value = nxlog.create_dataset(name='value', data=[property.value])

    if value and property.units:
        value.attrs['units'] = property.units

    # TODO should get times but current examples don't have that
    try:
        times = property.times
        if len(times) > 0:
            # convert to float in seconds
            epoch = times[0].toISO8601String()
            times = (times - times[0]) / np.timedelta64(1, 's')
            times = nxlog.create_dataset(name='time', data=times)
            times.attrs['offset'] = np.string_(epoch)
            times.attrs['units'] = 'second'
    except AttributeError:
        pass  # doesn't have times

    return nxlog


def _savespecialparameters(nxentry, wksp, **kwargs):
    nxcollection = _createnxgroup(nxentry, 'derived_parameters',
                                  'NXcollection')
    runObj = mtd[str(wksp)].run()

    # TODO add more of the "important" parameters
    names = ['beam_center_x', 'beam_center_y',
             'qmax', 'qmin', 'qstep', ]
    for name in names:
        if name in runObj:  # they are optional
            _savenxlog(nxcollection, runObj[name])


def savereductionlog(wksp1d, wksp2d=None, filename=None, **kwargs):
    r'''Save the reduction log TODO this needs more blah-blah describing the format

    Parameters
    ----------
    wksp1d: Workspace2D
        Workspace containing only one spectru (the I(q) curve)
    wksp2d: Workspace2D
        Not currently used
    filename: string
        The output filename to write to
    kwargs: dict
        dictionary of optional items:
        - ``python`` the script used to create everything
        - `pythonfile`` the name of the file containing the python script.
          Will be read into ``python`` argument if not already supplied
        - ``starttime`` when the original script was started
    '''
    if not filename:
        raise RuntimeError('Cannot write to file "{}"'.format(filename))
    if not wksp1d or str(wksp1d) not in mtd:
        raise RuntimeError('Cannot write out 1d workspace "{}"'.format(wksp1d))
    # TODO more checks?

    # save the 1d dataset
    SaveNexusProcessed(InputWorkspace=str(wksp1d), Filename=filename,
                       PreserveEvents=False)

    # TODO save the 2d data

    # re-open the file to append other information
    with h5py.File(filename, 'a') as handle:
        entry = _createnxgroup(handle, 'reduction_information', 'NXentry')

        # read the contents of the script
        _savepythonscript(entry, **kwargs)
        _savereductionparams(entry, **kwargs)

        # timestamp of when it happened - default to now
        starttime = kwargs.get('starttime', datetime.now().isoformat())
        entry.create_dataset(name='start_time', data=[np.string_(starttime)])

        # computer it was on
        hostname = environ.get('HOSTNAME', '')
        hostname = kwargs.get('hostname', hostname)  # parameter wins
        if hostname:
            entry.create_dataset(name='hostname', data=[np.string_(hostname)])

        # software involved
        _savenxprocess(entry, 'mantid', mantid_version)
        _savenxprocess(entry, 'sangria', sangria_version)

        # user information
        user = environ.get('USER', '')
        user = kwargs.get('user', user)
        if user:
            nxuser = _createnxgroup(entry, 'user', 'NXuser')
            nxuser.create_dataset(name='facility_user_id',
                                  data=[np.string_(user)])

            username = environ.get('USERNAME', '')
            username = kwargs.get('username', username)
            if username:
                nxuser.create_dataset(name='name',
                                      data=[np.string_(username)])

        _savespecialparameters(entry, wksp1d)

        # TODO   - add the logs of stdout and stderr
