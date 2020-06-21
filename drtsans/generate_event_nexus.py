# Create Event NeXus file
import drtsans
from drtsans.h5_buffer import HDFNode


class BankNode(drtsans.h5_buffer.GroupNode):
    """Node for bank entry such as /entry/bank12

    """
    def __init__(self):
        """

        """
        super(BankNode, self).__init__()

    def set_events(self):
        """

        Returns
        -------

        """
        pass


class InstrumentNode(drtsans.h5_buffer.GroupNode):
    """
    Node for instrument entry (i.e., /entry/instrument)
    """
    def __init__(self):
        """

        """
        super(InstrumentNode, self).__init__()


class DasLogNode(drtsans.h5_buffer.GroupNode):
    """
    Node for one specific DAS log such as /entry/DASlogs/sample_detector_distance
    """
    def __init__(self, log_name, log_times, log_values):
        """

        Parameters
        ----------
        log_name: str
            full path log name as /entry/DASlogs/{log_name}
        log_times
        log_values
        """
        super(DasLogNode, self).__init__(name=log_name)
        self._log_times = log_times
        self._log_values = log_values


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
