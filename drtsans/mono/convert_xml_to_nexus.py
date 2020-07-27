# This module contains the class and static method to convert legacy
# CG2's data file in XML format to event Nexus
import os
import numpy as np
from drtsans.mono.spice_xml_parser import SpiceXMLParser
from drtsans.files.event_nexus_rw import DasLog, EventNeXusWriter, TofHistogram
import h5py
from mantid.simpleapi import LoadHFIRSANS
from mantid.simpleapi import mtd   # logger


class EventNexusConverter(object):
    """
    Class to provide service to convert to event NeXus from various input
    """
    def __init__(self, beam_line, instrument_name):
        """

        Parameters
        ----------
        beam_line: str
            beam line such as CG2, CG3
        instrument_name: str
            instrument name such as CG2, CG3
        """
        # beam line name
        self._beam_line = beam_line
        self._instrument_name = instrument_name

        # instrument XML IDF content
        self._idf_content = None

        # counts
        self._detector_counts = None   # 1D array ranging from PID 0 (aka workspace index 0)
        self._monitor_counts = None

        # sample logs
        self._das_logs = dict()

        # run time
        self._run_start = None
        self._run_stop = None

        # run number
        self._run_number = None

    def generate_event_nexus(self, target_nexus, num_banks):
        """Generate event NeXus properly

        Parameters
        ----------
        target_nexus: str
            name of the output Nexus file
        num_banks: int
            CG2 = 48, CG3 = ???

        Returns
        -------

        """
        # Set constants
        pulse_duration = 0.1  # second
        tof_min = 1000.
        tof_max = 20000.

        # Generate event nexus writer
        event_nexus_writer = EventNeXusWriter(beam_line=self._beam_line, instrument_name=self._instrument_name)

        # set instrument
        event_nexus_writer.set_instrument_info(num_banks,  self._idf_content)

        # set counts
        for bank_id in range(1, num_banks + 1):
            # create TofHistogram instance
            start_pid, end_pid = self.get_pid_range(bank_id)
            pix_ids = np.arange(start_pid, end_pid + 1)
            counts = self._detector_counts[start_pid:end_pid + 1]
            counts = counts.astype('int64')
            histogram = TofHistogram(pix_ids, counts, pulse_duration, tof_min, tof_max)

            # set to writer
            event_nexus_writer.set_bank_histogram(bank_id, histogram)

        # set meta
        for das_log in self._das_logs.values():
            event_nexus_writer.set_meta_data(das_log)

        # set fake run number
        event_nexus_writer.set_run_number(self._run_number)

        # Write file
        event_nexus_writer.generate_event_nexus(target_nexus, self._run_start, self._run_stop, self._monitor_counts)

    def load_sans_xml(self, xml_file_name, prefix=''):
        """Load data and meta data from legacy SANS XML data file

        Parameters
        ----------
        xml_file_name: str
            name of SANS XML file
        prefix: str
            prefix for output workspace name

        Returns
        -------

        """
        # Load meta data and convert to NeXus format
        spice_log_dict, pt_number = self.retrieve_meta_data(xml_file_name)
        self._das_logs = self.convert_log_units(spice_log_dict)

        # output workspace name
        sans_ws_name = os.path.basename(xml_file_name).split('.xml')[0]
        sans_ws_name = f'{prefix}{sans_ws_name}'

        # load
        LoadHFIRSANS(Filename=xml_file_name,
                     OutputWorkspace=sans_ws_name)

        # get counts and reshape to (N, )
        sans_ws = mtd[sans_ws_name]
        counts = sans_ws.extractY().transpose().reshape((sans_ws.getNumberHistograms(),))

        self._detector_counts = counts[2:]
        monitor_counts = int(counts[0])
        self._monitor_counts = monitor_counts

        # get run start time: force to a new time
        self._run_start = sans_ws.run().getProperty('run_start').value
        self._run_stop = sans_ws.run().getProperty('end_time').value

        self._run_number = pt_number

    def load_idf(self, template_nexus_file):
        """Load IDF content from a template NeXus file

        Parameters
        ----------
        template_nexus_file

        Returns
        -------

        """
        # Import source
        source_nexus_h5 = h5py.File(template_nexus_file, 'r')
        # IDF in XML
        self._idf_content = source_nexus_h5['entry']['instrument']['instrument_xml']['data'][0]
        # Close
        source_nexus_h5.close()

        return

    @staticmethod
    def retrieve_meta_data(spice_file_name):
        """Retrieve meta from workspace

        Parameters
        ----------
        spice_file_name: str
            full path of SPICE data file in XML format

        Returns
        -------
        tuple
            (1) dict with key: Nexus das log name, value: (log value, log unit)  (2) int for pt number

        """
        # 'attenuator': 'attenuator_pos',
        # TODO FIXME This is very instrument-dependent!
        # <dcal pos="5.00000" units="mm" description="Detector Calibration Bar (mm)" type="FLOAT32">5.000005</dcal>
        # FIXME - shall source-distance written as  source-aperture-sample-aperture
        das_spice_log_map = {'CG2:CS:SampleToSi': ('sample_to_flange', 'mm', float),   # same
                             'sample_detector_distance': ('sdd', 'm', float),  # same
                             'wavelength': ('lambda', 'angstroms', float),  # angstroms -> A
                             'wavelength_spread': ('dlambda', 'fraction', float),  # fraction -> None
                             'source_aperture_diameter': ('source_aperture_size', 'mm', float),  # same
                             'sample_aperture_diameter': ('sample_aperture_size', 'mm', float),  # same
                             'detector_trans_Readback': ('detector_trans', 'mm', float),  # same
                             'source_distance': ('source_distance', 'm', float),  # same
                             'beamtrap_diameter': ('beamtrap_diameter', 'mm', float),   # not there
                             'dcal_Readback': ('dcal', 'mm', float)  # same unit mm
                             }

        # Load SPICE file
        spice_reader = SpiceXMLParser(spice_file_name)

        das_log_values = dict()
        for nexus_log_name, spice_tuple in das_spice_log_map.items():
            # read value from XML node
            spice_log_name, default_unit, data_type = spice_tuple
            value, unit = spice_reader.get_node_value(spice_log_name, data_type)

            # set default
            if unit is None:
                unit = default_unit
            else:
                # check unit
                if unit != default_unit:
                    raise RuntimeError(f'SPICE log {spice_log_name} has unit {unit} different from '
                                       f'expected {default_unit}')

            das_log_values[nexus_log_name] = value, unit

        # Attenuator is special
        das_log_values['attenuator'] = spice_reader.read_attenuator()

        # Get pt number
        pt_number, unit = spice_reader.get_node_value('Scan_Point_Number', int)

        spice_reader.close()

        return das_log_values, pt_number

    @staticmethod
    def convert_log_units(spice_log_dict):
        """Convert log unit from SPICE log to Nexus log

        Parameters
        ----------
        spice_log_dict:  ~dict
            key: Nexus das log name, value: (log value, log unit)

        Returns
        -------
         ~dict
            key: DAS log name, value: DasLog

        """
        nexus_log_dict = dict()

        for nexus_das_log_name in spice_log_dict:
            # get value
            log_value, log_unit = spice_log_dict[nexus_das_log_name]
            # use the name of the NeXus das log value unit
            if nexus_das_log_name == 'wavelength':
                log_unit = 'A'
            elif nexus_das_log_name == 'wavelength_spread':
                log_unit = None
            # form das log
            nexus_das_log = DasLog(nexus_das_log_name, np.array([0.]), np.array([log_value]), log_unit, None)
            # add
            nexus_log_dict[nexus_das_log_name] = nexus_das_log

        return nexus_log_dict

    @staticmethod
    def get_pid_range(bank_id):
        """Set GPSANS bank and pixel ID relation

        Parameters
        ----------
        bank_id: int
            bank ID from 1 to 48

        Returns
        -------
        tuple
            start PID, end PID (assuming PID are consecutive in a bank and end PID is inclusive)

        """
        # calculate starting PID
        if bank_id <= 24:
            # from 1 to 24: front panel
            start_pid = (bank_id - 1) * 2 * 1024
        else:
            # from 25 to 48: back panel
            start_pid = ((bank_id - 25) * 2 + 1) * 1024

        # calculate end PID
        end_pid = start_pid + 1023

        return start_pid, end_pid
