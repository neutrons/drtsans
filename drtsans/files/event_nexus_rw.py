import numpy as np
from collections import namedtuple

__all__ = ['TofHistogram', 'NexusEvents', 'EventNeXusWriter', 'generate_events_from_histogram',
           'convert_events_to_histogram']

# Specify parameter
# Histogram converted from TOF events
TofHistogram = namedtuple('TofHistogram', ['pixel_ids', 'counts', 'pulse_duration', 'tof_min', 'tof_max'])

# TOF events generated from histogram
NexusEvents = namedtuple('NexusEvents', ['event_id', 'event_index', 'event_time_offset', 'event_time_zero'])


class EventNeXusWriter(object):
    """
    Write an Event NeXus file
    """
    def __init__(self):
        """ Initialization
        """
        # Bank of events
        self._banks_dict = dict()

        # Meta data
        self._meta_data_dict = dict()

        # Run start time
        self._run_start = None

    def set_counts(self, bank_id, counts, detector_ids):
        self._banks_dict[bank_id] = counts, detector_ids

    def set_meta_data(self, meta_name, value, unit):
        self._meta_data_dict[meta_name] = value, unit

    def set_run_start_time(self, run_start_time):
        self._run_start = run_start_time


def generate_events_from_histogram(bank_histogram, tof_resolution=0.1):
    """Convert histogram (counts on detector pixels) to 'fake' events

    Parameters
    ----------
    bank_histogram: TofHistogram
        Histogram for a single bank
    tof_resolution: float
        resolution of TOF

    Returns
    -------
    NexusEvents
        Generated TOF events in NeXus format

    """
    # get total counts
    total_counts = bank_histogram.counts.sum()

    # Create event_id
    event_id_array = np.ndarray(shape=(total_counts,), dtype='uint32')
    # repeat each pixel for its 'counts' times to simulate the number of events
    start_index = 0
    for pid_index, pixel_id in enumerate(bank_histogram.pixel_ids):
        # get counts
        stop_index = start_index + bank_histogram.counts[pid_index]
        # set repetition
        event_id_array[start_index:stop_index] = pid_index
        # promote to next round
        start_index = stop_index

    # Get pulse related parameters
    num_events_per_pulse = int((bank_histogram.tof_max - bank_histogram.tof_min) / tof_resolution)
    num_pulses = total_counts // num_events_per_pulse  # This value is just a whole number. It could +1

    # event_time_offset, event_index
    single_pulse_tof = np.arange(num_events_per_pulse) * 0.1 + bank_histogram.tof_min
    event_time_offset_array = np.repeat(single_pulse_tof, num_pulses)
    # event indexes: number of events of each pulse: same value for the first N pulses completely filled
    event_index_array = np.zeros((num_pulses, ), dtype='uint64') + num_events_per_pulse
    # event_time_zero: range as [0, 1, ....] * pulse duration. such as [0, 0.16, 0.33, ...]
    event_time_zero_array = np.arange(num_pulses) * bank_histogram.pulse_duration

    # for the rest of events in the REAL last pulse: partially filled
    last_pulse_event_number = total_counts - num_pulses * num_events_per_pulse
    if last_pulse_event_number > 0:
        num_pulses += 1
        # add the incomplete TOF
        event_time_offset_array = np.concatenate((event_time_offset_array,
                                                  single_pulse_tof[0:last_pulse_event_number]))
        # add one more pulse
        if len(event_time_zero_array) > 0:
            prev_last_pulse_time = event_time_zero_array[-1]
        else:
            prev_last_pulse_time = 0
        last_pulse_time = prev_last_pulse_time + bank_histogram.pulse_duration
        event_time_zero_array = np.concatenate((event_time_zero_array, np.array([last_pulse_time])))
        # append number of events in last pulse
        event_index_array = np.concatenate((event_index_array, np.array([last_pulse_event_number])))

    # construct output
    faked_nexus_events = NexusEvents(event_id_array, event_index_array,
                                     event_time_offset_array, event_time_zero_array)

    return faked_nexus_events


def convert_events_to_histogram(bank_entry):
    """Convert events information in bank to histogram

    Parameters
    ----------
    bank_entry:

    Returns
    -------
    TofHistogram
        Histogram converted from TOF bank information

    """

    # Calculate counts on each pixel:w
    bank_event_ids = bank_entry['event_id']

    # Count each value's appearance in bank event ids array
    # result is from 0 to largest pixel ID in the event_id array
    per_pixel_counts = np.bincount(bank_event_ids)
    # ignore those pixels without any counts
    pixel_id_array = np.where(per_pixel_counts > 0)[0]
    pixel_counts_array = per_pixel_counts[per_pixel_counts > 0]

    # Get pulse time duration information
    pulse_start_times = bank_entry['event_time_zero'][()]
    pulse_duration = (pulse_start_times[1:] - pulse_start_times[:-1]).mean()

    # Get reasonable TOF information
    tof_array = bank_entry['event_time_offset'][()]

    # Define namedtuple to record histogram from TOF events
    histogram = TofHistogram(pixel_id_array, pixel_counts_array, pulse_duration,
                             tof_array.min(), tof_array.max())

    return histogram
