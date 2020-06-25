import numpy as np
from collections import namedtuple


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


def convert_histogram_to_events(det_id_array, det_counts_array, pulse_duration,
                                min_tof=2000, max_tof=1000, tof_resolution=0.1):
    """Convert histogram (counts on detector pixels) to 'fake' events

    Parameters
    ----------
    det_id_array
    det_counts_array
    pulse_duration: float
        pulse period duration in unit of second
    min_tof: float
        minimum TOF value in unit of microsecond
    max_tof: float
        maximum TOF value in unit of microsecond

    Returns
    -------
    ~tuple
        event_id (array), event_index (array), pulse_time_offset (array), event_time_zero, total_counts

    """
    # get total counts

    return


def convert_to_histogram_bank(bank_entry):
    """Convert events information in bank to histogram

    Parameters
    ----------
    bank_entry

    Returns
    -------
    ~namedtuple
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
    TofHistogram = namedtuple('Histogram', ['pixel_ids', 'counts', 'pulse_duration', 'tof_min', 'tof_max'])
    histogram = TofHistogram(pixel_id_array, pixel_counts_array, pulse_duration,
                             tof_array.min(), tof_array.max())

    return histogram
