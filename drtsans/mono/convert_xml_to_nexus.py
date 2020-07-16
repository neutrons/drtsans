# This module contains the class and static method to convert legacy
# CG2's data file in XML format to event Nexus
import os
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
        self._detector_counts = None
        self._monitor_counts = None

        # sample logs
        self._das_logs = dict()

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
        self._monitor_counts = counts[:2]

    @staticmethod
    def retrieve_meta_data(workspace):
        """Retrieve meta from workspace

        Parameters
        ----------
        workspace: ~mantid.Workspace2D
            workspace to retrieve meta data from

        Returns
        -------

        """
        das_spice_log_map = {'CG2:CS:SampleToSi': 'sample_to_flange',
                             'sample_detector_distance': 'sdd',
                             'wavelength': 'lambda',
                             'wavelength_spread': 'dlambda',
                             'source_aperture_diameter': 'source_aperture_size',
                             'sample_aperture_diameter': 'sample_aperture_size',
                             'detector_trans': 'detector_trans',
                             'attenuator': 'attenuator_pos',
                             'source_distance': 'source_distance',
                             'beamtrap_diameter': 'beamtrap_diameter'
                             }

        das_log_dict = dict()

        for nexus_log_name, spice_log_name in das_spice_log_map.items():
            try:
                p = workspace.run().getProperty(spice_log_name)
            except RuntimeError:
                print(f'{spice_log_name} cannot be found')
            else:
                das_log_dict[nexus_log_name] = p.times, p.value,  p.units
