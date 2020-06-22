# Create Event NeXus file
import numpy as np
import drtsans
from drtsans.h5_buffer import HDFNode, DataSetNode


class BankNode(drtsans.h5_buffer.GroupNode):
    """Node for bank entry such as /entry/bank12

    """
    def __init__(self):
        """

        """
        super(BankNode, self).__init__()

        # add NX_class
        self.add_attributes({'NX_class': 'NXevent_data'})

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
        super(InstrumentNode, self).__init__(name='/entry/instrument')

        # add the NeXus class attributes
        self.add_attributes({'NX_class': 'NXinstrument'})

    def set_instrument_info(self, target_station_number, beam_line, name):
        """

        Parameters
        ----------
        target_station_number: int
            target station number.  1 is used for HFIR
        beam_line: str
            CG2, CG3
        name: str
            CG2, CG3

        Returns
        -------
        None

        """
        str_single_array = np.ndarray(shape=(1,), dtype=np.dtype.str)

        # target station node
        target_station_node = DataSetNode(name=f'{self.name}/target_station_number')
        target_station_node.set_value(np.array(target_station_number))

        # beam line
        beam_line_node = DataSetNode(name=f'{self.name}/beamline')
        str_single_array[0] = beam_line
        beam_line_node.set_value(str_single_array)

        # beam line name
        name_node = DataSetNode(name=f'{self.name}/name')
        str_single_array[0] = name
        name_node.set_value(str_single_array)


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

        self.add_attributes({'NX_class': 'NXlog'})


class DasLogsCollectionNode(drtsans.h5_buffer.GroupNode):
    """
    Node for '/entry/DASlogs'
    """
    def __init__(self):
        """
        Initialization
        """
        super(DasLogsCollectionNode, self).__init__(name='/entry/DASlogs')
        self.add_attributes({'NX_class': 'NXcollection'})


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
