import pytest
import os
import h5py
from drtsans.files.event_nexus_rw import EventNeXusWriter
from drtsans.files.event_nexus_rw import generate_events_from_histogram
from drtsans.files.event_nexus_rw import convert_events_to_histogram


def test_write_event_nexus():
    assert EventNeXusWriter


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
    bank9_histogram = convert_events_to_histogram(bank9_entry)
    #  pixel_id_array, pixel_counts_array, pulse_duration, tof_array.min(), tof_array.max()

    # close  file
    nexus_h5.close()

    # verify
    assert bank9_histogram.counts.shape == (992, )
    # pixel ID 17000 (ws index 17000 too) at index 597 has 224 counts
    assert bank9_histogram.pixel_ids[597] == 17000
    assert bank9_histogram.counts[597] == 224
    assert bank9_histogram.pixel_ids.min() >= 16384 and bank9_histogram.pixel_ids.max() < 17408
    assert bank9_histogram.pulse_duration == pytest.approx(0.01666667, 1.E-4)
    assert bank9_histogram.tof_min >= 0.0
    assert bank9_histogram.tof_max == pytest.approx(16666.2, 0.1)


def test_convert_histogram_to_events(reference_dir):
    """

    Returns
    -------

    """
    # Create a TofHistogram from bank9
    nexus_name = os.path.join(reference_dir.new.gpsans, 'CG2_9166.nxs.h5')
    nexus_h5 = h5py.File(nexus_name, 'r')
    # test with bank 9
    bank9_entry = nexus_h5['/entry/bank9_events']
    bank9_histogram = convert_events_to_histogram(bank9_entry)
    # close  file
    nexus_h5.close()

    # generate events
    nexus_events = generate_events_from_histogram(bank9_histogram, 0.1)

    # Verification
    # check number of events
    assert nexus_events.event_id.shpae[0] == nexus_events.event_index.sum()
    assert nexus_events.event_id.shape == nexus_events.event_time_offset.shape
    assert nexus_events.event_index.shape == nexus_events.event_time_zero.shape
    assert nexus_events.event_time_offset.min() == bank9_histogram.tof_min
    assert nexus_events.event_time_offset.max() <= bank9_histogram.tof_max
    assert nexus_events.event_id.shape[0] == bank9_entry['total_counts'][0]


if __name__ == '__main__':
    pytest.main(__file__)
