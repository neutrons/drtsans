import pytest
import os
import h5py
from drtsans.files.event_nexus_rw import EventNeXusWriter
from drtsans.files.event_nexus_rw import convert_histogram_to_events
from drtsans.files.event_nexus_rw import convert_to_histogram_bank


def test_imports():
    assert EventNeXusWriter


def test_convert_histogram_to_events():
    """

    Returns
    -------

    """
    assert convert_histogram_to_events


def test_convert_to_histogram(reference_dir):
    """Test method to convert a single bank's TOF events to histogram

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Parse NeXus file manually for the values nodes
    nexus_name = os.path.join(reference_dir.new.gpsans, 'CG2_9166.nxs.h5')
    nexus_h5 = h5py.File(nexus_name, 'r')

    # test with bank 9
    bank9_entry = nexus_h5['/entry/bank9_events']
    bank9_histogram = convert_to_histogram_bank(bank9_entry)
    #  pixel_id_array, pixel_counts_array, pulse_duration, tof_array.min(), tof_array.max()

    # close  file
    nexus_h5.close()

    # verify
    assert bank9_histogram.counts.shape == (982, )
    assert bank9_histogram.pixel_ids.min() >= 16384 and bank9_histogram.pixel_ids.max() < 17408
    assert bank9_histogram.pulse_duration == pytest.approx(0.01666667, 1.E-6)
    assert bank9_histogram.min_tof >= 0.0
    assert bank9_histogram.max_tof == pytest.approx(16666.6, 0.1)


if __name__ == '__main__':
    pytest.main(__file__)
