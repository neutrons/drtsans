import h5py
import numpy as np
import pytest
import os

from drtsans import savereductionlog
from drtsans.iq import determine_1d_log_bins
from tests.unit.new.drtsans.i_of_q_binning_tests_data import generate_test_data, get_gold_1d_log_bins
from drtsans.dataobjects import IQmod
from drtsans.dataobjects import IQazimuthal
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


def _create_iq():
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

    # Test the high level method
    test_iq = IQmod(intensities, sigmas, scalar_q_array, scalar_dq_array)

    return test_iq


def _create_iqxqy():
    # Get data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)
    test_iqxqy = IQazimuthal(intensity=intensities,
                             error=sigmas,
                             qx=qx_array,
                             qy=qy_array,
                             delta_qx=dqx_array,
                             delta_qy=dqy_array)
    return test_iqxqy


def _create_tmp_log_filename():
    tmp_log_filename = NamedTemporaryFile(prefix="logfile", suffix='.nxs.h5').name
    tmp_log_filename = os.path.abspath(tmp_log_filename)
    if os.path.exists(tmp_log_filename):
        os.remove(tmp_log_filename)
    return tmp_log_filename


def _checkNXData(nxentry, name):
    nxdata = _getGroup(nxentry, name, 'NXdata')

    print(_strValue(nxdata, 'data'))

    assert False

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


def test_writing_metadata():

    pythonscript = "this is my python script"
    reductionparams = "reduction parameter 1: value1"
    starttime = '1993-03-18T21:00:00'
    username = 'Neymar'
    user = 'Cavani'

    test_iq = _create_iq()
    tmp_log_filename = _create_tmp_log_filename()
    savereductionlog(tmp_log_filename, iq=test_iq,
                     python=pythonscript,
                     starttime=starttime,
                     user=user,
                     username=username,
                     reductionparams=reductionparams)

    assert os.path.exists(tmp_log_filename), 'log file {} does not exist'.format(tmp_log_filename)

    # with h5py.File(tmp_log_filename, 'r') as handle:
    #     _checkProcessingEntry(handle,
    #                           pythonscript=pythonscript,
    #                           starttime=starttime,
    #                           reductionparams=reductionparams,
    #                           user=user,
    #                           username=username)

def _test_data(tested_data=[], ref_data=[], abs=None):
    for _tested, _ref in zip(tested_data, ref_data):
        if abs is None:
            assert _tested == _ref
        else:
            _tested == pytest.approx(_ref, abs=abs)


def test_writing_iq():
    test_iq = _create_iq()
    tmp_log_filename = _create_tmp_log_filename()
    savereductionlog(tmp_log_filename, iq=test_iq)

    assert os.path.exists(tmp_log_filename), 'log file {} does not exist'.format(tmp_log_filename)

    with h5py.File(tmp_log_filename, 'r') as handle:
        iq_nxdata = _getGroup(handle, 'I(Q)', 'NXdata')
        data = iq_nxdata['I'][:]
        _test_data(tested_data=data,
                   ref_data=np.array([93, 60]))

        data = iq_nxdata['Idev'][:]
        _test_data(tested_data=data,
                   ref_data=np.array([9.64365076, 7.74596669]),
                   abs=1e-7)

        data = iq_nxdata['Q'][:]
        _test_data(tested_data=data,
                   ref_data=np.array([0.0078897, 0.0059338]),
                   abs=1e-7)

        data = iq_nxdata['Qdev'][:]
        _test_data(tested_data=data,
                   ref_data=np.array([0.011912, 0.11912]),
                   abs=1e-6)

        # assert False



def test_writing_iqxqy():
    test_iqxqy = _create_iqxqy()
    tmp_log_filename = _create_tmp_log_filename()
    savereductionlog(tmp_log_filename, iqxqy=test_iqxqy)

    assert os.path.exists(tmp_log_filename), 'log file {} does not exist'.format(tmp_log_filename)

    with h5py.File(tmp_log_filename, 'r') as handle:
        pass



def test_writing_iq_and_iqxqy():
    test_iq = _create_iq()
    test_iqxqy = _create_iqxqy()
    tmp_log_filename = _create_tmp_log_filename()
    savereductionlog(tmp_log_filename, iq=test_iq, iqxqy=test_iqxqy)

    assert os.path.exists(tmp_log_filename), 'log file {} does not exist'.format(tmp_log_filename)

    with h5py.File(tmp_log_filename, 'r') as handle:
        pass


def test_no_data_passed():
    with pytest.raises(RuntimeError):
        savereductionlog()


if __name__ == '__main__':
    pytest.main([__file__])
