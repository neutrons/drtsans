# Unit test for converting SPICE XML data file to event NeXus
import pytest
import h5py
import os
import numpy as np
from drtsans.mono.gpsans.cg2_spice_to_nexus import CG2EventNexusConvert
from drtsans.mono.biosans.cg3_spice_to_nexus import CG3EventNexusConvert


def test_cg2_pid_range(reference_dir):
    """Test PID range

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Load test event NeXus file
    test_nexus = os.path.join(reference_dir.new.gpsans, 'CG2_9177.nxs.h5')
    nexus_h5 = h5py.File(test_nexus, 'r')

    # Check each bank
    for bank_id in range(1, 48 + 1):  # 48 banks
        pids = nexus_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        min_pid = np.min(pids)
        max_pid = np.max(pids)
        start_pid, end_pid = CG2EventNexusConvert().get_pid_range(bank_id)

        assert start_pid <= min_pid <= max_pid <= end_pid, f'CG2 Bank {bank_id} PID is out of range'

    # close file
    nexus_h5.close()

    # Check out of range bank ID
    with pytest.raises(RuntimeError):
        CG2EventNexusConvert().get_pid_range(0)
    with pytest.raises(RuntimeError):
        CG2EventNexusConvert().get_pid_range(49)


def test_cg3__pid_range(reference_dir):
    """Test PID range

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Load test event NeXus file
    test_nexus = os.path.join(reference_dir.new.biosans, 'CG3_5705.nxs.h5')
    nexus_h5 = h5py.File(test_nexus, 'r')

    # Check each bank
    for bank_id in range(1, 88 + 1):  # 88 banks
        pids = nexus_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        min_pid = np.min(pids)
        max_pid = np.max(pids)
        start_pid, end_pid = CG3EventNexusConvert().get_pid_range(bank_id)

        assert start_pid <= min_pid <= max_pid <= end_pid, f'CG3 Bank {bank_id} PID (H5 range ' \
                                                           f'{min_pid} - {max_pid}) is out of range ' \
                                                           f'of calculated PID range {start_pid} - ' \
                                                           f'{end_pid}'

    # close file
    nexus_h5.close()

    # Check out of range bank ID
    with pytest.raises(RuntimeError):
        CG3EventNexusConvert().get_pid_range(0)
    with pytest.raises(RuntimeError):
        CG3EventNexusConvert().get_pid_range(89)


if __name__ == '__main__':
    pytest.main(__file__)
