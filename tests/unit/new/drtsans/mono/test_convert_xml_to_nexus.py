# Unit test for converting SPICE XML data file to event NeXus
import pytest
import h5py
import os
import numpy as np
from drtsans.mono.gpsans.cg2_spice_to_nexus import CG2EventNexusConvert


def test_pid_range(reference_dir):
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

        assert start_pid <= min_pid <= max_pid <= end_pid, f'Bank {bank_id} PID is out of range'

    # close file
    nexus_h5.close()


if __name__ == '__main__':
    pytest.main(__file__)
