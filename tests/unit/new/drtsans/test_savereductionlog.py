import h5py
from mantid.simpleapi import mtd, CompareWorkspaces, Load, LoadNexusProcessed
import numpy as np
from drtsans import savereductionlog
import os
import pytest
from tempfile import gettempdir, NamedTemporaryFile


def _strValue(group, name):
    '''Get a value from a SDS'''
    assert name in group, 'Did not find "{}" in "{}"'.format(name, group.name)
    return group[name].value[0].decode('utf-8')


def _strAttr(data, name):
    '''Get an attribute as a string'''
    assert name in data.attrs, \
        '"{}" does not have attribute "{}"'.format(data.name, name)
    value = data.attrs[name]
    try:
        return value.decode('utf-8')
    except AttributeError:
        return value


def _getGroup(parent, name, klass):
    '''Utility function to get a group in the hdf5 file. The default error messages are lacking.'''
    assert name in parent, \
        '{} does not contain a group named "{}" '.format(parent.name, name)
    child = parent[name]
    assert _strAttr(child, 'NX_class') == klass, \
        '{} is not of type "{}" '.format(name, klass)
    return child


def _check1d(handle, wksp_name):
    '''Utility function for verifying the 1d data (and attributes) are correct'''
    wksp = mtd[wksp_name]
    dims = (wksp.getNumberHistograms(), wksp.blocksize())

    # TODO should there be a better name for the entry?
    entry = _getGroup(handle, 'mantid_workspace_1', 'NXentry')

    assert _strValue(entry, 'workspace_name') == wksp_name

    nxdata = entry['workspace']

    axis1 = nxdata['axis1']
    assert axis1.size == dims[1]+1  # for a histogram
    assert _strAttr(axis1, 'units') == 'MomentumTransfer'
    assert np.all(axis1.value == wksp.readX(0))

    axis2 = nxdata['axis2']
    assert axis2.size == dims[0]
    assert _strAttr(axis2, 'units') == 'spectraNumber'
    assert axis2.value == 1.

    values = nxdata['values']
    assert np.all(values.shape == dims)
    assert _strAttr(values, 'units') == 'Counts'
    assert _strAttr(values, 'axes') == 'axis2,axis1'
    assert values.attrs['signal'] == 1
    assert np.all(values.value == wksp.readY(0))

    errors = nxdata['errors']
    assert np.all(errors.value == wksp.readE(0))


def _checkNXNote(nxentry, name, mimetype, file_name, data):
    '''Utility function for verifying that the NXnote has the
    appropriate information'''
    # TODO question: should these be there when file_name
    # and data are both empty?
    nxnote = _getGroup(nxentry, name, 'NXnote')

    assert _strValue(nxnote, 'type') == mimetype
    assert _strValue(nxnote, 'file_name') == file_name
    assert _strValue(nxnote, 'data') == data


def _checkNXprocess(entry, program):
    '''Utility function for verifying that the NXprocess has the
    appropriate information'''
    nxprocess = _getGroup(entry, program, 'NXprocess')
    assert _strValue(nxprocess, 'program') == program
    assert _strValue(nxprocess, 'version')  # having one is enough


def _checkNXcollection(nxentry, name, param_names):
    '''Utility function for verifying that the NXcollection has the
    appropriate information'''
    nxcollection = _getGroup(nxentry, name, 'NXcollection')
    for param in param_names:
        assert param in nxcollection, \
            'Could not find "{}" in "{}"'.format(param, nxcollection.name)
        nxlog = nxcollection[param]
        assert 'value' in nxlog, '{} did not have a value'.format(param)

        if 'time' in nxlog:
            times = nxlog['time']
            assert _strAttr(times, 'units') == 'second'
            assert _strAttr(times, 'offset')  # that there is one is enough


def _checkProcessingEntry(handle, **kwargs):
    '''Utility function for verifying that the processing NXentry has the
    appropriate information'''
    entry = _getGroup(handle, 'reduction_information', 'NXentry')

    if 'starttime' in kwargs:
        assert 'start_time' in entry
        assert _strValue(entry, 'start_time') == kwargs['starttime']
    assert 'hostname' in entry

    nxuser = _getGroup(entry, 'user', 'NXuser')
    assert 'facility_user_id' in nxuser
    assert _strValue(nxuser, 'name') == kwargs['username']

    _checkNXNote(entry, 'reduction_script', 'text/x-python',
                 kwargs.get('pythonfile', ''),
                 kwargs.get('pythonscript', ''))
    _checkNXNote(entry, 'reduction_parameters',
                 'application/json',
                 '', kwargs.get('reductionparams', ''))
    param_names = ['beam_center_x', 'beam_center_y']
    _checkNXcollection(entry, 'derived_parameters', param_names)
    _checkNXprocess(entry, 'mantid')
    _checkNXprocess(entry, 'drtsans')


def _checkWorkspaces(filename, orig, entry):
    '''Utility function for verifying that the workspace saved is the
    same as the one that is in the file'''
    if not orig:
        print('nothing to check against')
        return

    reloaded = orig + '_reload'
    LoadNexusProcessed(Filename=filename, OutputWorkspace=reloaded,
                       EntryNumber=entry)
    result, msg = CompareWorkspaces(Workspace1=orig,
                                    Workspace2=reloaded)
    assert result, msg
    if reloaded in mtd:
        mtd.remove(reloaded)


@pytest.mark.parametrize(
    'filename1d', ['test_save_output/EQSANS_68200_iq.nxs', '']
)
@pytest.mark.parametrize(
    'filename_other', [(),
                       ('test_save_output/EQSANS_68200_iq.nxs', ''),
                       ('', 'test_save_output/EQSANS_68200_iq.nxs'),
                       ('test_save_output/EQSANS_68200_iq.nxs',
                        'test_save_output/EQSANS_68200_iq.nxs')]
)
def test_saving(reference_dir, filename1d, filename_other):
    wksp1d = ''
    wksp_other = []

    # setup inputs
    if filename1d:
        filename1d = os.path.join(reference_dir.new.eqsans,
                                  filename1d)
        wksp1d = 'test_save_wksp1d'
        Load(Filename=filename1d, OutputWorkspace=wksp1d)

    for i, filename in enumerate(filename_other):
        if filename:
            filename = os.path.join(reference_dir.new.eqsans,
                                    filename)
            wksp = 'test_save_wksp_{}'.format(i)
            if filename == filename1d:
                wksp = wksp1d  # just reuse the workspace
            else:
                Load(Filename=filename, OutputWorkspace=wksp)
            wksp_other.append(wksp)
        else:
            wksp_other.append('')

    tmpfile = NamedTemporaryFile(prefix=wksp1d, suffix='.nxs.h5').name
    tmpfile = os.path.abspath(tmpfile)
    if os.path.exists(tmpfile):
        os.unlink(tmpfile)  # remove it if it already exists

    # dummy arguments to check against - they should be found in the file
    pythonscript = 'blah blah blah'
    reductionparams = ''
    starttime = '1992-01-19T00:00:01Z'
    username = 'Jimmy Neutron'
    user = 'neutron'

    # run the function - use same workspace for both
    if wksp1d:
        savereductionlog(tmpfile, wksp1d, *wksp_other,
                         python=pythonscript, starttime=starttime,
                         user=user, username=username)

        # validation
        assert os.path.exists(tmpfile), \
            'Output file "{}" does not exist'.format(tmpfile)

        # go into the file to check things
        with h5py.File(tmpfile, 'r') as handle:
            _check1d(handle, wksp1d)
            _checkProcessingEntry(handle, pythonscript=pythonscript,
                                  reductionparams=reductionparams,
                                  starttime=starttime, username=username)

        # use mantid to check the workspaces
        _checkWorkspaces(tmpfile, wksp1d, 1)
        for i, wksp in enumerate(wksp_other):
            if wksp:
                _checkWorkspaces(tmpfile, wksp, i + 1)
    else:  # not supplying 1d workspace should always fail
        with pytest.raises(RuntimeError):
            savereductionlog(tmpfile, wksp1d, *wksp_other,
                             python=pythonscript, starttime=starttime,
                             user=user, username=username)
        assert not os.path.exists(tmpfile), \
            'Output file "{}" should not exist'.format(tmpfile)

    # cleanup
    if os.path.exists(tmpfile):
        os.unlink(tmpfile)
    wksp_other.append(wksp1d)
    for wksp in wksp_other:
        if wksp and wksp in mtd:
            mtd.remove(wksp)


def test_no_arguments():
    with pytest.raises(RuntimeError):
        savereductionlog()


def test_empty_filename():
    with pytest.raises(RuntimeError):
        savereductionlog(filename='', wksp=None)


def test_nonexistant_1d_wksp():
    tmpfile = os.path.join(gettempdir(), 'test_nonexistant_1d_wksp.nxs')
    with pytest.raises(RuntimeError):
        savereductionlog(filename=tmpfile, wksp=None)


if __name__ == '__main__':
    pytest.main([__file__])
