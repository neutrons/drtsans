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
    '''Write an NXnote inside an entry
    http://download.nexusformat.org/doc/html/classes/base_classes/NXnote.html

    nxentry: HDF handle
        NXentry to put the NXnote into
    name: str
        Name of the NXnote group
    mimetype: str
        Mimetype of the data
    file_name: str
        Name of the file that the note was taken from
    data: str
        The contents of the note
    '''
    nxnote = _createnxgroup(nxentry, name, 'NXnote')

    nxnote.create_dataset(name='file_name', data=[np.string_(file_name)])
    nxnote.create_dataset(name='type', data=[np.string_(mimetype)])
    nxnote.create_dataset(name='data', data=[np.string_(data)])

    return nxnote


def _savepythonscript(nxentry, pythonfile, pythonscript):
    '''Write the python script as a NXnote

    nxentry: HDF handle
        NXentry to put the NXnote into
    pythonfile: str
        Filename that was supplied to reduction
    pythonscript: str
        The python script itself
    '''
    # read the contents of the script
    if not pythonscript:
        with open(pythonfile, 'r') as script:
            pythonscript = script.read()

    return _savenxnote(nxentry, 'reduction_script', 'text/x-python',
                       pythonfile, pythonscript)


def _savereductionparams(nxentry, parameters):
    '''Save the all of the reduction parameters as an NXnote

    nxentry: HDF handle
        NXentry to put the NXnote into
    parameters: dict or str
        The parameters supplied to the reduction script. This will be converted
        to a json string if it isn't one already.
    '''
    if not isinstance(parameters, str):
        parameters = json.dumps(parameters)
    return _savenxnote(nxentry, 'reduction_parameters',
                       'application/json', file_name='',
                       data=parameters)


def _savenxprocess(nxentry, program, version):
    '''Create a NXprocess for the specified program

    Parameters
    ----------
    nxentry: HDF handle
        NXentry to put the NXprocess into
    program: str
        name of the program
    version: str
        program's version string
    '''
    nxprocess = _createnxgroup(nxentry, program, 'NXprocess')
    nxprocess.create_dataset(name='program', data=[np.string_(program)])
    nxprocess.create_dataset(name='version', data=[np.string_(version)])


def _savenxlog(nxcollection, property):
    '''Create a NXlog from the supplied property

    Parameters
    ----------
    nxcollection: HDF handle
        NXcollection that the NXlog should be added to
    property: PropertyWithValue
        log item to inspect and write its values to disk
    '''
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


def _savespecialparameters(nxentry, wksp):
    '''Save the special parameters from the workspace

    Currently this only saves the derived parameters. It should be expanded to
    save some of the special input parameters as well.

    Parameters
    ----------
    nxentry: HDF handle
        Entry group to put information in
    wksp: Workspace2D
        Workspace to get the information from
    '''
    nxcollection = _createnxgroup(nxentry, 'derived_parameters',
                                  'NXcollection')
    runObj = mtd[str(wksp)].run()

    # TODO add more of the "important" derived parameters
    names = ['beam_center_x', 'beam_center_y',
             'qmax', 'qmin', 'qstep', ]
    for name in names:
        if name in runObj:  # they are optional
            _savenxlog(nxcollection, runObj[name])

    # TODO look in the workspace history for some of the special input
    #      parameters


def savereductionlog(wksp, filename, *args, **kwargs):
    r'''Save the reduction log

    There are three ``NXentry``. The first is for the 1d reduced data, second
    is for the 2d reduced data, and the third is for the extra information
    about how the data was processed.

    Parameters
    ----------
    wksp: Workspace2D
        Workspace containing only one spectru (the I(q) curve). This is the
        only workspace that metadata is taken from.
    filename: string
        The output filename to write to
    args: Workspaces
        Other workspaces to save to the file. This can be the Workspace
        containing I(Qx, Qy), or multiple 1d Workspaces for additional frames,
        or a mixture
    kwargs: dict
        dictionary of optional items:
        - ``python`` the script used to create everything
        - `pythonfile`` the name of the file containing the python script.
          Will be read into ``python`` argument if not already supplied
        - ``reductionparams`` is the parameters supplied to the reduction
          script
        - ``starttime`` when the original script was started
        - ``hostname`` name of the computer used. If not provided, will be
          gotten from the system environment ``HOSTNAME``
        - ``user`` user-id of who reduced the data (as in xcamms). If not
          provided will be gotten from the system environment ``USER``
        - ``username`` username of who reduced the data (as in actual name).
          If not provided will be gotten from the system environment
          ``USERNAME``
    '''
    if not filename:
        raise RuntimeError('Cannot write to file "{}"'.format(filename))
    if not wksp or str(wksp) not in mtd:
        raise RuntimeError('Cannot write out 1d workspace "{}"'.format(wksp))
    # TODO more checks?

    # save the 2 workspaces
    SaveNexusProcessed(InputWorkspace=str(wksp), Filename=filename,
                       PreserveEvents=False)
    for optional_wksp in args:
        if optional_wksp:
            SaveNexusProcessed(InputWorkspace=str(optional_wksp),
                               Filename=filename,
                               PreserveEvents=False, Append=True)

    # re-open the file to append other information
    with h5py.File(filename, 'a') as handle:
        entry = _createnxgroup(handle, 'reduction_information', 'NXentry')

        # read the contents of the script
        _savepythonscript(entry, pythonfile=kwargs.get('pythonfile', ''),
                          pythonscript=kwargs.get('python', ''))
        _savereductionparams(entry, parameters=kwargs.get('reductionparams',
                                                          ''))

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

        _savespecialparameters(entry, wksp)

        # TODO   - add the logs of stdout and stderr
