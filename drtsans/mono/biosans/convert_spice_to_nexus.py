from drtsans.mono.convert_xml_to_nexus import EventNeXusWriter
from drtsans.files.event_nexus_rw import parse_event_nexus


DAS_LOGs = ['CG3:CS:SampleToSi', 'sample_detector_distance', 'wavelength', 'wavelength_spread',
            'source_aperture_diameter', 'sample_aperture_diameter',
            'detector_trans_Readback', 'ww_rot_Readback',
            'source_aperture_sample_aperture_distance']


def generate_event_nexus(source_nexus, target_nexus, das_log_list=DAS_LOGs):
    """Generate event NeXus properly from a source Nexus file

    This method will be migrated to drtsans.mono.biaosans

    Parameters
    ----------
    source_nexus: str
        source nexus file
    target_nexus: str
        target nexus file
    das_log_list: ~list
        list of DAS logs

    Returns
    -------

    """
    cg3_num_banks = 88

    # Import essential experimental data from source event nexus file
    nexus_contents = parse_event_nexus(source_nexus, 88, das_log_list)
    # Generate event nexus writer
    event_nexus_writer = EventNeXusWriter(beam_line='CG3', instrument_name='CG3')

    # set instrument: 88 banks (2 detectors)
    event_nexus_writer.set_instrument_info(cg3_num_banks,  nexus_contents[0])

    # set counts: 88 banks (2 detectors)
    for bank_id in range(1, cg3_num_banks + 1):
        event_nexus_writer.set_bank_histogram(bank_id, nexus_contents[1][bank_id])

    # set meta
    for das_log in nexus_contents[5].values():
        event_nexus_writer.set_meta_data(das_log)

    # time
    start_time = nexus_contents[3]
    end_time = nexus_contents[4]

    # Write file
    event_nexus_writer.generate_event_nexus(target_nexus, start_time, end_time, nexus_contents[2])
