import h5py
from mantid.simpleapi import mtd, CompareWorkspaces, Load, LoadNexusProcessed
import numpy as np
from ornl.sans import savereductionlog
import os
import pytest
from tempfile import gettempdir


def strValue(group, name):
    assert name in group, 'Did not find "{}" in "{}"'.format(name, group.name)
    return group[name].value[0].decode('utf-8')


def strAttr(data, name):
    assert name in data.attrs, \
        '"{}" does not have attribute "{}"'.format(data.name, name)
    value = data.attrs[name]
    try:
        return value.decode('utf-8')
    except AttributeError:
        return value


def getGroup(parent, name, klass):
    assert name in parent, \
        '{} does not contain a group named "{}" '.format(parent.name, name)
    child = parent[name]
    assert strAttr(child, 'NX_class') == klass, \
        '{} is not of type "{}" '.format(name, klass)
    return child


def check1d(handle, wksp_name):
    wksp = mtd[wksp_name]
    dims = (wksp.getNumberHistograms(), wksp.blocksize())

    # TODO should there be a better name for the entry?
    entry = getGroup(handle, 'mantid_workspace_1', 'NXentry')

    assert strValue(entry, 'workspace_name') == wksp_name

    nxdata = entry['workspace']

    axis1 = nxdata['axis1']
    assert axis1.size == dims[1]+1  # for a histogram
    assert strAttr(axis1, 'units') == 'MomentumTransfer'
    assert np.all(axis1.value == wksp.readX(0))

    axis2 = nxdata['axis2']
    assert axis2.size == dims[0]
    assert strAttr(axis2, 'units') == 'spectraNumber'
    assert axis2.value == 1.

    values = nxdata['values']
    assert np.all(values.shape == dims)
    assert strAttr(values, 'units') == 'Counts'
    assert strAttr(values, 'axes') == 'axis2,axis1'
    assert values.attrs['signal'] == 1
    assert np.all(values.value == wksp.readY(0))

    errors = nxdata['errors']
    assert np.all(errors.value == wksp.readE(0))


def checkNXNote(nxentry, name, mimetype, file_name, data):
    # TODO question: should these be there when file_name
    # and data are both empty?
    nxnote = getGroup(nxentry, name, 'NXnote')

    assert strValue(nxnote, 'type') == mimetype
    assert strValue(nxnote, 'file_name') == file_name
    assert strValue(nxnote, 'data') == data


def checkNXprocess(entry, program):
    nxprocess = getGroup(entry, program, 'NXprocess')
    assert strValue(nxprocess, 'program') == program
    assert strValue(nxprocess, 'version')  # having one is enough


def checkNXcollection(nxentry, name, param_names):
    nxcollection = getGroup(nxentry, name, 'NXcollection')
    for param in param_names:
        assert param in nxcollection, \
            'Could not find "{}" in "{}"'.format(param, nxcollection.name)
        nxlog = nxcollection[param]
        assert 'value' in nxlog, '{} did not have a value'.format(param)

        if 'time' in nxlog:
            times = nxlog['time']
            assert strAttr(times, 'units') == 'second'
            assert strAttr(times, 'offset')  # that there is one is enough


def checkProcessingEntry(handle, **kwargs):
    entry = getGroup(handle, 'reduction_information', 'NXentry')

    if 'starttime' in kwargs:
        assert 'start_time' in entry
        assert strValue(entry, 'start_time') == kwargs['starttime']
    assert 'hostname' in entry

    nxuser = getGroup(entry, 'user', 'NXuser')
    assert 'facility_user_id' in nxuser
    assert strValue(nxuser, 'name') == kwargs['username']

    checkNXNote(entry, 'reduction_script', 'text/x-python',
                kwargs.get('pythonfile', ''),
                kwargs.get('pythonscript', ''))
    checkNXNote(entry, 'reduction_parameters',
                'application/json',
                '', kwargs.get('reductionparams', ''))
    param_names = ['beam_center_x', 'beam_center_y']
    checkNXcollection(entry, 'derived_parameters', param_names)
    checkNXprocess(entry, 'mantid')
    checkNXprocess(entry, 'sangria')


def checkWorkspaces(filename, orig, entry):
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
    'filename2d', ['test_save_output/EQSANS_68200_iq.nxs', '']
)
def test_saving(refd, filename1d, filename2d):
    wksp1d = ''
    wksp2d = ''

    # setup inputs
    tmpfile = os.path.join(gettempdir(), wksp1d + '.nxs.h5')
    tmpfile = os.path.abspath(tmpfile)
    if os.path.exists(tmpfile):
        os.unlink(tmpfile)  # remove it if it already exists

    if filename1d:
        filename1d = os.path.join(refd.new.eqsans,
                                  filename1d)
        wksp1d = 'test_save_wksp1d'
        Load(Filename=filename1d, OutputWorkspace=wksp1d)

    if filename2d:
        filename2d = os.path.join(refd.new.eqsans,
                                  filename2d)
        wksp2d = 'test_save_wksp2d'
        if filename2d == filename1d:
            wksp2d = wksp1d  # just reuse the workspace
        else:
            Load(Filename=filename2d, OutputWorkspace=wksp2d)

    pythonscript = 'blah blah blah'
    reductionparams = ''
    starttime = '1992-01-19T00:00:01Z'
    username = 'Jimmy Neutron'
    user = 'neutron'

    # run the function - use same workspace for both
    if wksp1d:
        savereductionlog(wksp1d, wksp2d, filename=tmpfile, python=pythonscript,
                         starttime=starttime, user=user, username=username)

        # validation
        assert os.path.exists(tmpfile), \
            'Output file "{}" does not exist'.format(tmpfile)

        # go into the file to check things
        with h5py.File(tmpfile, 'r') as handle:
            check1d(handle, wksp1d)
            checkProcessingEntry(handle, pythonscript=pythonscript,
                                 reductionparams=reductionparams,
                                 starttime=starttime, username=username)

        # use mantid to check the workspaces
        checkWorkspaces(tmpfile, wksp1d, 1)
        checkWorkspaces(tmpfile, wksp2d, 2)
    else:  # not supplying 1d workspace should always fail
        with pytest.raises(RuntimeError):
            savereductionlog(wksp1d, wksp2d, filename=tmpfile,
                             python=pythonscript, starttime=starttime,
                             user=user, username=username)
        assert not os.path.exists(tmpfile), \
            'Output file "{}" should not exist'.format(tmpfile)

    # cleanup
    if os.path.exists(tmpfile):
        os.unlink(tmpfile)
    for wksp in [wksp1d, wksp2d]:
        if wksp and wksp in mtd:
            mtd.remove(wksp)


def test_empty_filename():
    with pytest.raises(RuntimeError):
        savereductionlog(None, filename='')


def test_nonexistant_1d_wksp():
    tmpfile = os.path.join(gettempdir(), 'test_nonexistant_1d_wksp.nxs')
    with pytest.raises(RuntimeError):
        savereductionlog(None, filename=tmpfile)


if __name__ == '__main__':
    pytest.main()
