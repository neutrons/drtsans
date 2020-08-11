# Convert GPSANS SPICE file to event NeXus
import os
from drtsans.mono.gpsans.cg2_spice_to_nexus import CG2EventNexusConvert

# Template event nexus file for instrument geometry and etc
TEMPLATE_EVENT_NEXUS = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/gpsans/CG2_9177.nxs.h5'


def convert_spice_to_nexus(ipts_number, exp_number, scan_number, pt_number):
    """ Convert one SPICE to NeXus

    Parameters
    ----------
    ipts_number: int
        IPTS
    exp_number: int
        experiment number
    scan_number: int
        scan
    pt_number: int
        pt

    Returns
    -------

    """
    # Build the path to the data
    spice_dir = f'/HFIR/CG2/IPTS-{ipts_number}/exp{exp_number}/Datafiles'
    assert os.path.exists(spice_dir), f'SPICE data directory {spice_dir} cannot be found'
    spice_data_file = os.path.join(spice_dir,
                                   f'CG2_exp{exp_number}_scan{scan_number:04}_{pt_number:04}.xml')
    assert os.path.exists(spice_data_file), f'SPICE file {spice_data_file} cannot be located'

    # Template Nexus file
    template_nexus_file = TEMPLATE_EVENT_NEXUS
    assert os.path.exists(template_nexus_file), f'Template NeXus file {template_nexus_file} cannot be located'

    # Specify the output directory
    output_dir = f'/HFIR/CG2/IPTS-{ipts_number}/shared/spice_nexus'
    if not os.path.exists(output_dir):
        raise RuntimeError(f'Output NeXus directory {output_dir} does not exist.'
                           f'Create directory {output_dir} and grand access to all IPTS users')

    # output file name
    out_nexus_file = os.path.join(output_dir, f'CG2_{exp_number}{scan_number:04}{pt_number:04}.nxs.h5')
    out_nexus_file = os.path.join(output_dir, out_nexus_file)

    # Load meta data and convert to NeXus format
    das_log_map = {'CG2:CS:SampleToSi': ('sample_to_flange', 'mm'),  # same
                   'sample_detector_distance': ('sdd', 'm'),  # same
                   'wavelength': ('lambda', 'angstroms'),  # angstroms -> A
                   'wavelength_spread': ('dlambda', 'fraction'),  # fraction -> None
                   'source_aperture_diameter': ('source_aperture_size', 'mm'),  # same
                   'sample_aperture_diameter': ('sample_aperture_size', 'mm'),  # same
                   'detector_trans_Readback': ('detector_trans', 'mm'),  # same
                   'source_distance': ('source_distance', 'm'),  # same. source-aperture-sample-aperture
                   'beamtrap_diameter': ('beamtrap_diameter', 'mm'),  # not there
                   'attenuator': ('attenuator_pos', 'mm')  # special
                   }

    # init converter
    converter = CG2EventNexusConvert()
    # load instrument definition (IDF)
    converter.load_idf(template_nexus_file)
    # load SPICE (xml file)
    converter.load_sans_xml(spice_data_file, das_log_map)
    # generate event nexus
    converter.generate_event_nexus(out_nexus_file)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


# set
ipts = 828
exp = 280
scan = 4
pt = 1

convert_spice_to_nexus(ipts, exp, scan, pt)