# Create Event NeXus file
import numpy as np
import drtsans
from drtsans.h5_buffer import DataSetNode, GroupNode


class BankNode(drtsans.h5_buffer.GroupNode):
    """Node for bank entry such as /entry/bank12

    """
    def __init__(self):
        """

        """
        super(BankNode, self).__init__()

        # add NX_class
        self.add_attributes({'NX_class': b'NXevent_data'})

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
        self.add_attributes({'NX_class': b'NXinstrument'})

    def set_instrument_info(self, target_station_number, beam_line, name, short_name):
        """

        Parameters
        ----------
        target_station_number: int
            target station number.  1 is used for HFIR
        beam_line: Bytes
            CG2, CG3
        name: Bytes
            CG2, CG3

        Returns
        -------
        None

        """
        # target station node
        target_station_node = DataSetNode(name=f'{self.name}/target_station_number')
        target_station_node.set_value(np.array(target_station_number))
        self.set_child(target_station_node)

        # beam line
        beam_line_node = DataSetNode(name=f'{self.name}/beamline')
        beam_line_node.set_1d_string([beam_line])
        self.set_child(beam_line_node)

        # beam line name
        name_node = DataSetNode(name=f'{self.name}/name')
        name_node.set_1d_string([name])
        self.set_child(name_node)
        name_node.add_attributes({'short_name': short_name})

    def set_idf(self, idf_str, idf_type, description):
        """Set instrument xml

        Parameters
        ----------
        idf_str: Bytes
            IDF XML string
        idf_type: Bytes
            IDF type
        description: Bytes
            Description

        Returns
        -------

        """
        # Create the instrument_xml node
        xml_node_name = f'{self.name}/instrument_xml'
        xml_node = GroupNode(name=xml_node_name)
        xml_node.add_attributes({'NX_class': b'NXnote'})
        self.set_child(xml_node)

        # add data node
        data_node = DataSetNode(name=f'{xml_node_name}/data')
        data_node.set_1d_string([idf_str])
        xml_node.set_child(data_node)

        # add description
        des_node = DataSetNode(name=f'{xml_node_name}/description')
        des_node.set_1d_string([description])
        xml_node.set_child(des_node)

        # add type
        type_node = DataSetNode(name=f'{xml_node_name}/type')
        type_node.set_1d_string([idf_type])
        xml_node.set_child(type_node)


class DasLogNode(drtsans.h5_buffer.GroupNode):
    """
    Node for one specific DAS log such as /entry/DASlogs/sample_detector_distance
    """
    def __init__(self, log_name, log_times, start_time, log_values, log_unit):
        """DAS log node for specific

        Parameters
        ----------
        log_name: str
            full path log name as /entry/DASlogs/{log_name}
        log_times: numpy.ndarray
            relative sample log time
        start_time: str
            ISO standard time for run start
        log_values: numpy.ndarray
            sample log values
        log_unit: Byes
            log unit
        """
        super(DasLogNode, self).__init__(name=log_name)
        self._log_times = log_times
        self._run_start = start_time
        self._log_values = log_values
        self._log_unit = log_unit

        self.add_attributes({'NX_class': b'NXlog'})

        self._set_time_value()

    def _set_time_value(self):
        """Set time and value including
        - average_value
        - average_value_error
        - maximum_value
        - minimum_value
        - time
        - value

        Returns
        -------

        """

        return

    def set_device_info(self, device_id, device_name, target):
        """Set node for device related information

        Parameters
        ----------
        device_id
        device_name
        target

        Returns
        -------

        """
        # Create Device ID node
        for node_name, info_value in [('device_id', device_id),
                                      ('device_name', device_name),
                                      ('target', target)]:
            child_node = DataSetNode(name=self._create_child_name(node_name))
            child_node.set_value(np.array(info_value))
            self._children.append(child_node)

        return


class DasLogsCollectionNode(drtsans.h5_buffer.GroupNode):
    """
    Node for '/entry/DASlogs'
    """
    def __init__(self):
        """
        Initialization
        """
        super(DasLogsCollectionNode, self).__init__(name='/entry/DASlogs')
        self.add_attributes({'NX_class': b'NXcollection'})


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
