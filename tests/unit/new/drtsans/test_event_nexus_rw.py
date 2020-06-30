import pytest
import numpy as np
import os
import h5py
from drtsans.files.event_nexus_rw import EventNeXusWriter
from drtsans.files.event_nexus_rw import generate_events_from_histogram, generate_monitor_events_from_count
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
    total_counts = bank9_entry['total_counts'][0]
    # close  file
    nexus_h5.close()

    # generate events
    nexus_events = generate_events_from_histogram(bank9_histogram, 0.1)

    # Verification
    # event index only contain the starting event index of each pulse. Its aggregated value is useless
    assert nexus_events.event_id.shape[0] >= nexus_events.event_index.sum()
    assert nexus_events.event_id.shape == nexus_events.event_time_offset.shape
    assert nexus_events.event_index.shape == nexus_events.event_time_zero.shape
    assert nexus_events.event_time_offset.min() == bank9_histogram.tof_min
    assert nexus_events.event_time_offset.max() <= bank9_histogram.tof_max
    # check number of events:
    assert nexus_events.event_id.shape[0] == total_counts


def test_convert_monitor_counts_to_events():
    """Test convert monitor counts to events

    Returns
    -------

    """
    # Generate pulse time: 0.16 second for 250 pulses
    event_time_zero_array = np.arange(250).astype('float') * 0.16

    tof_min = 0.
    tof_max = 10000

    # Case 1: Counts can be spread evenly
    num_counts = 25000

    # Execute
    monitor_events_even = generate_monitor_events_from_count(num_counts, event_time_zero_array, tof_min, tof_max)

    # Verify
    assert monitor_events_even.event_time_offset.shape == (num_counts, )
    assert monitor_events_even.event_index.shape == (250, )
    assert monitor_events_even.event_index[-1] == num_counts - 100

    # Case 2: Counts cannot be spread unevenly
    num_counts = 25013

    # Execute
    monitor_events_uneven = generate_monitor_events_from_count(num_counts, event_time_zero_array, tof_min, tof_max)

    # Verify
    assert monitor_events_uneven.event_time_offset.shape == (num_counts, )
    assert monitor_events_uneven.event_index.shape == (250, )
    assert monitor_events_uneven.event_index[-1] == num_counts - 101
    assert monitor_events_uneven.event_index[-12] == num_counts - 101 * 12


if __name__ == '__main__':
    pytest.main(__file__)
