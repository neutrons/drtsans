# import h5py
# from mantid.simpleapi import mtd, CompareWorkspaces, Load, LoadNexusProcessed
import numpy as np
from drtsans import savereductionlog
import pytest
import os
from drtsans.iq import determine_1d_log_bins
from tests.unit.new.drtsans.i_of_q_binning_tests_data import generate_test_data, get_gold_1d_log_bins
from drtsans.dataobjects import IQmod
from tempfile import NamedTemporaryFile


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


# def _checkWorkspaces(filename, orig, entry):
#     '''Utility function for verifying that the workspace saved is the
#     same as the one that is in the file'''
#     if not orig:
#         print('nothing to check against')
#         return
#
#     reloaded = orig + '_reload'
#     LoadNexusProcessed(Filename=filename, OutputWorkspace=reloaded,
#                        EntryNumber=entry)
#     result, msg = CompareWorkspaces(Workspace1=orig,
#                                     Workspace2=reloaded)
#     assert result, msg
#     if reloaded in mtd:
#         mtd.remove(reloaded)


def test_writing_iq():
    # Define Q range from tab '1D_bin_log_no_sub_no_wt' in r4
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    num_steps_per_10 = 10  # 10 steps per decade

    # Verify bin edges and bin center
    log_bins = determine_1d_log_bins(q_min, q_max, num_steps_per_10)
    gold_edges, gold_centers = get_gold_1d_log_bins()

    np.testing.assert_allclose(log_bins.edges, gold_edges, rtol=5.E-4)
    np.testing.assert_allclose(log_bins.centers, gold_centers, rtol=5.E-4)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    # binned_iq = _do_1d_no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
    #                                      log_bins.centers, log_bins.edges)

    # Test the high level method
    test_iq = IQmod(intensities, sigmas, scalar_q_array, scalar_dq_array)
    # binned_iq = bin_intensity_into_q1d(test_iq, log_bins, BinningMethod.NOWEIGHT)

    tmp_log_filename = NamedTemporaryFile(prefix="logfile", suffix='.nxs.h5').name
    tmp_log_filename = os.path.abspath(tmp_log_filename)
    if os.path.exists(tmp_log_filename):
        os.remove(tmp_log_filename)
    savereductionlog(tmp_log_filename, iq=test_iq)

    assert os.path.exists(tmp_log_filename)


def test_no_data_passed():
    with pytest.raises(RuntimeError):
        savereductionlog()


if __name__ == '__main__':
    pytest.main([__file__])
