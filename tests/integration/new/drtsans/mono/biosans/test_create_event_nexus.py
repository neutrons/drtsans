"""
Integration test to create event nexus file
"""
import pytest
import numpy as np
import os
from drtsans.load import load_events
import json
import h5py
from drtsans.files.hdf5_rw import GroupNode, DataSetNode
from drtsans.files.event_nexus_nodes import InstrumentNode, DasLogNode, BankNode, MonitorNode
from drtsans.files.event_nexus_rw import generate_events_from_histogram
from drtsans.files.event_nexus_rw import generate_monitor_events_from_count
from drtsans.files.event_nexus_rw import init_event_nexus, parse_event_nexus, EventNeXusWriter
# drtsans imports
from drtsans.mono.biosans import (load_all_files, reduce_single_configuration,
                                  reduction_parameters, validate_reduction_parameters)
from mantid.simpleapi import LoadEventNexus
from drtsans.files.hdf5_rw import FileNode


# FIXME - BioSANS special
# - dark current: das log 'duration'


def generate_event_nexus_prototype_x(source_nexus_file, prototype_dup_nexus):
    # Parse
    logs_white_list = ['CG3:CS:SampleToSi', 'sample_detector_distance',
                       'wavelength', 'wavelength_spread',
                       'source_aperture_diameter', 'sample_aperture_diameter',
                       'detector_trans_Readback']
    cg3_nexus = parse_event_nexus(source_nexus_file, 88, logs_white_list)

    # Load source
    source_h5 = h5py.File(source_nexus_file, 'r')
    source_root = FileNode()
    source_root.parse_h5_entry(source_h5)
    source_entry = source_root.get_child('/entry')

    # Create a new one
    duplicate_root = FileNode()
    duplicate_root.set_child(source_root.get_child('/entry'))
    duplicate_entry_node = duplicate_root.get_child('/entry')

    # Replace node instrument
    source_instrument = source_entry.get_child('/entry/instrument')
    source_xml = source_instrument.get_child('/entry/instrument/instrument_xml')
    xml_idf_content = source_xml.get_child('/entry/instrument/instrument_xml/data').value[0]

    new_instrument_node = InstrumentNode()
    new_instrument_node.set_instrument_info(1, 'CG3', 'CG3', 'CG3')
    new_instrument_node.set_idf(xml_idf_content,
                                idf_type='XML content of instrument IDF', description='text/xml')

    # The bank nodes
    bank_histograms = cg3_nexus[1]
    run_start_time = cg3_nexus[3]

    bank_node_dict = dict()
    for bank_id in range(1, 88 + 1):
        bank_node_name = f'/entry/bank{bank_id}_events'
        bank_node = duplicate_entry_node.get_child(bank_node_name)
        bank_node_dict[bank_id] = bank_node
        if bank_id not in [48, 53]:
            duplicate_entry_node.remove_child(bank_node_name)

    # Add back all the bank nodes
    max_pulse_time_array = None
    for bank_id in bank_node_dict:
        if bank_id in [48, 53]:
            nexus_events = generate_events_from_histogram(bank_histograms[bank_id], 10., verbose=True)
            print(f'Bank {bank_id}, Pixel IDs: {nexus_events.event_id.min()} to {nexus_events.event_id.max()}')
            continue
        # generate fake events from counts
        nexus_events = generate_events_from_histogram(bank_histograms[bank_id], 10.)
        # Create bank node for bank
        bank_node = BankNode(name=f'/entry/bank{bank_id}_events', bank_name=f'bank{bank_id}')
        bank_node.set_events(nexus_events.event_id, nexus_events.event_index,
                             nexus_events.event_time_offset, run_start_time,
                             nexus_events.event_time_zero)
        if max_pulse_time_array is None or nexus_events.event_time_zero.shape[0] > max_pulse_time_array.shape[0]:
            max_pulse_time_array = nexus_events.event_time_zero
        # set child
        duplicate_entry_node.set_child(bank_node)

    # Monitor counts
    duplicate_entry_node.remove_child('/entry/monitor1')
    tof_min = 0.
    tof_max = 10000.
    monitor_events = generate_monitor_events_from_count(cg3_nexus[2], max_pulse_time_array, tof_min, tof_max)
    target_monitor_node = MonitorNode('/entry/monitor1', 'monitor1')
    target_monitor_node.set_monitor_events(event_index_array=monitor_events.event_index,
                                           event_time_offset_array=monitor_events.event_time_offset,
                                           run_start_time=run_start_time,
                                           event_time_zero_array=max_pulse_time_array)
    duplicate_entry_node.set_child(target_monitor_node)

    # replace instrument node
    duplicate_entry_node.remove_child('/entry/instrument')
    duplicate_entry_node.set_child(new_instrument_node)

    # Delete nodes under entry
    black_list = ['title',
                  'total_counts',
                  'total_other_counts',
                  'total_pulses',
                  'total_uncounted_counts',
                  'user1',
                  'user2', 'notes',
                  'user3',
                  'user4',
                  'user5',
                  'entry_identifier', 'definition', 'bank_error_events', 'bank_unmapped_events',
                  'sample', 'experiment_title', 'experiment_identifier', 'run_number', 'proton_charge',
                  'raw_frames',
                  'Software'
                  ]
    for short_name in black_list:
        full_name = f'/entry/{short_name}'
        duplicate_entry_node.remove_child(full_name)

    # TODO  - Testing block
    # Sample logs
    das_logs_node = duplicate_entry_node.get_child('/entry/DASlogs')
    # get all the children names
    das_log_dict = dict()
    for child_log_node in das_logs_node.children:
        das_log_dict[child_log_node.name] = child_log_node

    # remove all children
    for child_log_name in das_log_dict:
        das_logs_node.remove_child(child_log_name)

    # add back some
    for child_log_name in das_log_dict:
        short_name = child_log_name.split('DASlogs/')[1]
        if short_name in logs_white_list:
            das_logs_node.set_child(das_log_dict[child_log_name])
    # END-TODO

    # write
    duplicate_root.write(prototype_dup_nexus)

    # Close
    source_h5.close()

    return


def test_copy_event_nexus(reference_dir):
    """Prototype test to find out why LoadEventNexusFiled

    LoadEventNexus-[Warning] Empty proton_charge sample log. You will not be able to filter by time.
    LoadEventNexus-[Error] Error in execution of algorithm LoadEventNexus:
    LoadEventNexus-[Error] Error finding workspace index; pixelID 49152 with offset 2 is out of range
        (length=49154)
    =================================================================================================

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Get the source file
    source_nexus_file = 'CG3_5709.nxs.h5'
    source_nexus_file = os.path.join(reference_dir.new.biosans, source_nexus_file)
    assert os.path.exists(source_nexus_file), f'Test data {source_nexus_file} does not exist'

    # Duplicate the source file to the temporary directory
    output_dir = '/tmp/prototype_cg3nexus'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    prototype_dup_nexus = os.path.join(output_dir, 'CG3_5709_prototype.nxs.h5')

    generate_event_nexus_prototype_x(source_nexus_file, prototype_dup_nexus)

    prototype_ws = load_events(prototype_dup_nexus, output_workspace='cg3_prototype', NumberOfBins=2)
    assert prototype_ws
    assert prototype_ws.getNumberHistograms() == 90112


def test_duplicate_event_nexus(reference_dir, cleanfile):
    """Test duplicating an HDF5/NeXus in 2 different approaches in order to verify EventNexusWriter

    Verification is to load both of the generated Event NeXus to do a comparison

    Test data: BIOSANS run 5709

    Returns
    -------

    """
    # Get the source file
    source_nexus_file = 'CG3_5709.nxs.h5'
    source_nexus_file = os.path.join(reference_dir.new.biosans, source_nexus_file)
    assert os.path.exists(source_nexus_file), f'Test data {source_nexus_file} does not exist'

    # Duplicate the source file to the temporary directory
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/dupcg3nexus'
    # cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    prototype_dup_nexus = os.path.join(output_dir, 'CG3_5709_prototype.nxs.h5')
    product_dup_nexus = os.path.join(output_dir, 'CG3_5709_product.nxs.h5')

    # Duplicate with both approach
    logs_white_list = ['CG3:CS:SampleToSi', 'sample_detector_distance',
                       'wavelength', 'wavelength_spread',
                       'source_aperture_diameter', 'sample_aperture_diameter',
                       'detector_trans_Readback']
    # generate_event_nexus_prototype(source_nexus_file, prototype_dup_nexus, logs_white_list)
    generate_event_nexus_prototype_x(source_nexus_file, prototype_dup_nexus)
    generate_event_nexus(source_nexus_file, product_dup_nexus, logs_white_list)

    # Load the duplicated
    prototype_ws = load_events(prototype_dup_nexus, output_workspace='cg3_prototype', NumberOfBins=2)

    # Load source file to workspace
    target_ws = load_events(product_dup_nexus, output_workspace='cg3_product', NumberOfBins=2)

    # Compare pixels' positions
    num_hist = prototype_ws.getNumberHistograms()
    for iws in range(0, num_hist, 100):
        source_det_i_pos = prototype_ws.getInstrument().getDetector(iws).getPos()
        target_det_i_pos = target_ws.getInstrument().getDetector(iws).getPos()
        np.testing.assert_allclose(source_det_i_pos, target_det_i_pos,
                                   err_msg=f'Mismatch is detected at Detector {iws}')
    # Check source position
    source_moderator_pos = prototype_ws.getInstrument().getSource().getPos()
    target_moderator_pos = target_ws.getInstrument().getSource().getPos()
    np.testing.assert_allclose(source_moderator_pos, target_moderator_pos,
                               err_msg=f'Mismatch is detected at neutron source position')

    # Compare counts on each pixel
    source_y = prototype_ws.extractY()
    target_y = target_ws.extractY()
    np.testing.assert_allclose(source_y, target_y)

    # Compare meta data
    if len(prototype_ws.getRun().getProperties()) != len(target_ws.getRun().getProperties()):
        print(f'Prototype:')
        for p in prototype_ws.getRun().getProperties():
            print(f'property: {p.name}')
        print(f'Product')
        for p in target_ws.getRun().getProperties():
            print(f'property: {p.name}')
        raise RuntimeError('Meta data mismatch')


def generate_event_nexus(source_nexus, target_nexus, das_log_list):
    """Generate event NeXus properly

    Parameters
    ----------
    source_nexus
    target_nexus

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


def generate_event_nexus_prototype(source_nexus, target_nexus, das_log_list):
    """Generate event NeXus using white list.

    This serves as the prototype to create event nexus from SANS histogram raw data

    White list
    Entry attributes: {'NX_class': b'NXentry'}

    White List Node: /entry/monitor1
    White List Node: /entry/proton_charge

    White List Node: /entry/duration
    White List Node: /entry/start_time
    White List Node: /entry/end_time

    White List Node: /entry/experiment_identifier
    White List Node: /entry/experiment_title
    White List Node: /entry/title
    White List Node: /entry/notes
    White List Node: /entry/raw_frames
    White List Node: /entry/run_number

    White List Node: /entry/total_counts
    White List Node: /entry/total_other_counts
    White List Node: /entry/total_pulses

    Parameters
    ----------
    source_nexus: str
        source event NeXus file name
    target_nexus

    Returns
    -------

    """
    # parse nexus information
    nexus_contents = parse_event_nexus(source_nexus, 88, das_log_list)

    # Create new nexus file structure
    target_nexus_root = init_event_nexus()

    target_entry_node = target_nexus_root.get_child('entry', is_short_name=True)

    # set instrument node
    set_instrument_node(nexus_contents[0], target_entry_node)

    # set DAS logs
    set_das_log_node(nexus_contents[5], nexus_contents[3], target_entry_node)

    # Add node on the white list
    entry_level_white_list = [
        ('/entry/start_time', nexus_contents[3]),
        ('/entry/end_time', nexus_contents[4])
    ]
    for child_node_name, child_value in entry_level_white_list:
        child_node = DataSetNode(child_node_name)
        child_node.set_string_value(child_value)
        target_entry_node.set_child(child_node)

    # set Bank 1 - 88 (2 detectors)
    max_pulse_time_array = None
    for bank_id in range(1, 88 + 1):
        bank_node_i = set_single_bank_node(nexus_contents[1][bank_id], target_entry_node, bank_id=bank_id,
                                           run_start_time=nexus_contents[3])
        event_time_zeros = bank_node_i.get_child('event_time_zero', is_short_name=True).value
        if max_pulse_time_array is None or event_time_zeros.shape[0] > max_pulse_time_array.shape[0]:
            max_pulse_time_array = event_time_zeros

    # Set monitor node
    set_monitor_node(nexus_contents[2],  nexus_contents[3], target_entry_node, max_pulse_time_array)

    # write
    target_nexus_root.write(target_nexus)


def set_monitor_node(monitor_counts, run_start_time, target_entry_node, event_time_zeros):
    """

    Parameters
    ----------
    monitor_counts: float, int
    target_entry_node
    event_time_zeros: ~numpy.ndarray
        event time zeros
    run_start_time: str, Bytes
        run start time

    Returns
    -------

    """
    # Generate a monitor node
    target_monitor_node = MonitorNode('/entry/monitor1', 'monitor1')

    tof_min = 0.
    tof_max = 10000.
    monitor_events = generate_monitor_events_from_count(monitor_counts, event_time_zeros, tof_min, tof_max)

    target_monitor_node.set_monitor_events(event_index_array=monitor_events.event_index,
                                           event_time_offset_array=monitor_events.event_time_offset,
                                           run_start_time=run_start_time,
                                           event_time_zero_array=event_time_zeros)

    target_entry_node.set_child(target_monitor_node)


def set_single_bank_node(bank_histogram, target_entry_node, bank_id, run_start_time):
    """Test writing bank 9 from histogram

    Parameters
    ----------
    bank_histogram: TofHistogram
        HDF5 file entry
    target_entry_node: GroupNode
        Target (output) group node for /entry/
    bank_id: int
        bank ID (from 1 to 88)
    run_start_time: str, Bytes
        run start time

    Returns
    -------
    BankNode
        newly generated bank node

    """
    # generate events
    nexus_events = generate_events_from_histogram(bank_histogram, 10.)

    try:
        run_start_time = np.string_(run_start_time).decode()
    except AttributeError:
        pass

    # Create bank node for bank
    bank_node = BankNode(name=f'/entry/bank{bank_id}_events', bank_name=f'bank{bank_id}')
    bank_node.set_events(nexus_events.event_id, nexus_events.event_index,
                         nexus_events.event_time_offset, run_start_time,
                         nexus_events.event_time_zero)

    # Link with parent
    target_entry_node.set_child(bank_node)

    return bank_node


def set_instrument_node(xml_idf, target_entry_node):
    """Set instrument node

    Parameters
    ----------
    xml_idf:  str
        IDF content
    target_entry_node

    Returns
    -------

    """
    # Create new instrument node
    instrument_node = InstrumentNode()
    target_entry_node.set_child(instrument_node)

    # Set values
    instrument_node.set_idf(xml_idf, idf_type=b'text/xml', description=b'XML contents of the instrument IDF')
    instrument_node.set_instrument_info(target_station_number=1, beam_line=b'CG3', name=b'CG3', short_name=b'CG3')


def set_das_log_node(das_log_dict, run_start_time, target_entry_node):
    """Set DAS log node in a mixed way

    Parameters
    ----------
    das_log_dict: dict
        das log dictionary containing DasLog objects
    run_start_time: str
        run start time
    target_entry_node: GroupNode
        target node

    Returns
    -------

    """
    target_logs_node = GroupNode('/entry/DASlogs')
    target_entry_node.set_child(target_logs_node)
    # add attribute
    target_logs_node.add_attributes({'NX_class': 'NXcollection'})

    for log_name in das_log_dict:
        set_single_log_node(target_logs_node, das_log_dict[log_name], run_start_time)


def set_single_log_node(log_collection_node, das_log, start_time):
    """

    Parameters
    ----------
    log_collection_node
    das_log: DasLog
    start_time: str

    Returns
    -------

    """
    # Set up a DAS log node
    das_log_node = DasLogNode(log_name=f'/entry/DASlogs/{das_log.name}',
                              log_times=das_log.times,
                              log_values=das_log.values,
                              start_time=start_time,
                              log_unit=das_log.unit)

    if das_log.device is not None:
        if das_log.device.target is None:
            device_target = None
        else:
            device_target = das_log.device.target
        das_log_node.set_device_info(device_id=das_log.device.id,
                                     device_name=das_log.device.name,
                                     target=device_target)

    # append to parent node
    log_collection_node.set_child(das_log_node)


def set_sdd_node(log_collection_node, source_h5):
    # Get times and value for /entry/DASlogs/sample_detector_distance
    ssd_entry = source_h5['entry']['DASlogs']['sample_detector_distance']
    ssd_times = ssd_entry['time'].value
    ssd_start_time = ssd_entry['time'].attrs['start']
    ssd_value = ssd_entry['value'].value
    ssd_value_unit = ssd_entry['value'].attrs['units']

    # Set up a DAS log node
    ssd_test_node = DasLogNode(log_name='/entry/DASlogs/sample_detector_distance',
                               log_times=ssd_times, log_values=ssd_value,
                               start_time=ssd_start_time, log_unit=ssd_value_unit)

    ssd_test_node.set_device_info(device_id=13, device_name=b'Mot-Galil3',
                                  target=b'/entry/DASlogs/CG3:CS:SampleToDetRBV')

    # append to parent node
    log_collection_node.set_child(ssd_test_node)


def verify_histogram(source_nexus, test_nexus):
    """Check whether two NeXus files can render out same result

    Parameters
    ----------
    source_nexus: str
        source/gold nexus file name
    test_nexus: str
        nexus file to test

    Returns
    -------

    """
    # Load NeXus file
    src_ws = LoadEventNexus(Filename=source_nexus, OutputWorkspace='gold', NumberOfBins=1)
    test_ws = LoadEventNexus(Filename=test_nexus, OutputWorkspace='test', NumberOfBins=1)

    # Compare counts
    error_message = ''
    for i in range(src_ws.getNumberHistograms()):
        if src_ws.readY(i)[0] != test_ws.readY(i)[0]:
            error_message += f'Workspace-index {i} / detector ID {src_ws.getDetector(i).getID()}/' \
                             f'{test_ws.getDetector(i).getID()}: Expected counts = {src_ws.readY(i)},' \
                             f'Actual counts = {test_ws.readY(i)}\n'

    # Nothing wrong!
    if error_message == '':
        return

    # Compare in HDF5 level
    source_h5 = h5py.File(source_nexus, 'r')
    target_h5 = h5py.File(test_nexus, 'r')

    # Check all the banks
    for bank_id in range(1, 89):
        source_event_ids = source_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        target_event_ids = target_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        if source_event_ids.shape != target_event_ids.shape:
            print(f'Bank {bank_id}  {source_event_ids.shape} vs {target_event_ids.shape}')
        for pid in [48139, 58367]:
            if source_event_ids.min() <= pid < source_event_ids.max():
                print(f'PID {pid} in source bank {bank_id}')
            if target_event_ids.min() <= pid < target_event_ids.max():
                print(f'PID {pid} in target bank {bank_id}')

    # close
    source_h5.close()
    target_h5.close()

    print(error_message)
    print(f'source: {source_nexus}')
    print(f'target: {test_nexus}')
    raise AssertionError(error_message)


def failed_test_reduction(reference_dir, cleanfile):
    """Test generate (partially copy) an event Nexus file by
    verifying reduction result between raw and generated event nexus file

    Testing is modified from mono.gpsans.test_overwrite_geometry_meta_data.test_no_overwrite()

    Returns
    -------

    """
    # Set up test
    json_str = generate_testing_json(os.path.join(reference_dir.new.biosans, 'overwrite_gold_04282020'), None, None)

    # TODO / FIXME - switch to tempfile later
    # output_dir = mkdtemp(prefix='meta_overwrite_bio_test1')
    output_dir = '/tmp/nexuscg3reduction/'
    # cleanfile(output_dir)

    # Run
    reduce_biosans_data(reference_dir.new.biosans, json_str, output_dir, prefix='BioMetaRaw')

    # Get result files
    sample_names = ['csmb_ecoli1h_n2']
    gold_path = os.path.join(reference_dir.new.biosans, 'overwrite_gold_04282020/test1/')

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, title='Raw (no overwriting)', prefix='test1')


def generate_testing_json(sens_nxs_dir, sample_to_si_window_distance, sample_to_detector_distance):
    """Generating testing JSON

    Parameters
    ----------
    sens_nxs_dir: str
        directory path to sensitivity files
    sample_to_si_window_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit millimeter
    sample_to_detector_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit meter

    Returns
    -------
    str
        JSON string
    """
    specs = {
        "iptsNumber": "23782",
        "sample": {"runNumber": "4822", "thickness": "0.1", "transmission": {"runNumber": "4822"}},
        "outputFileName": "CG3_4822",
        "background": {"runNumber": "4821", "transmission": {"runNumber": "4821"}},
        "beamCenter": {"runNumber": "1322"},
        "emptyTransmission": {"runNumber": "5705"},
        "configuration": {
            "outputDir": "/HFIR/CG3/shared/UserAcceptance/override/test1",
            "sampleApertureSize": "14",
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-12,244-256'}", "{'Bank':'21-24,45-48'}"],
            "darkMainFileName": "CG3_1383.nxs.h5",
            "darkWingFileName": "CG3_1383.nxs.h5",
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle486/sens_f4829m7p0_TDC_SAC.h5",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle486/sens_f4835w3p2_TDC_SAC.h5",
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": "0.0055e-8",
            "numMainQxQyBins": "100",
            "numWingQxQyBins": "100",
            "1DQbinType": "scalar",
            "QbinType": "log",
            "useLogQBinsEvenDecade": False,
            "LogQBinsPerDecadeMain": 20,
            "LogQBinsPerDecadeWing": 25,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "1",
            "QminMain": "0.003",
            "QminWing": "0.003",
            "overlapStitchQmin": "0.075",
            "overlapStitchQmax": "0.095",
            "useTimeSlice": False,
            "timeSliceInterval": "200",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'BIOSANS', validate=False)  # add defaults and defer validation
    reduction_config = reduction_input['configuration']  # a handy shortcut

    #  '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4829m7p0_TDC_SAC.h5'
    main_sens = os.path.join(sens_nxs_dir, 'sens_f4829m7p0_TDC_SAC.h5')
    reduction_config['sensitivityMainFileName'] = main_sens

    # '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4835w3p2_TDC_SAC.h5'
    wing_sens = os.path.join(sens_nxs_dir, 'sens_f4835w3p2_TDC_SAC.h5')
    reduction_config['sensitivityWingFileName'] = wing_sens

    if sample_to_si_window_distance is not None:
        reduction_config['sampleToSi'] = sample_to_si_window_distance

    if sample_to_detector_distance is not None:
        reduction_config['sampleDetectorDistance'] = sample_to_detector_distance

    return json.dumps(reduction_input)  # return a string representation


def reduce_biosans_data(nexus_dir, json_str, output_dir, prefix):
    """Reduce BIOSANS runs

    Parameters
    ----------
    nexus_dir: str
        path to NeXus files
    json_str: str
        configuration json
    output_dir: str
        output directory
    prefix: str
        prefix for output workspaces' names

    Returns
    -------

    """
    # Set up (testing) runs
    sample_names = ['csmb_ecoli1h_n2']
    sample = '5709'
    samples_tran = sample
    backgrounds = ['5715']
    backgrounds_trans = backgrounds

    # Replace duplicate
    source_sample_nexus = os.path.join(nexus_dir, f'CG3_{sample}.nxs.h5')
    os.path.exists(source_sample_nexus), f'Source sample NeXus {source_sample_nexus} does not exist'
    # TODO - need a better name for test sample nexus
    test_sample_nexus = os.path.join(output_dir, f'CG3_{sample}.nxs.h5')
    generate_event_nexus_prototype_x(source_sample_nexus, test_sample_nexus)
    verify_histogram(source_sample_nexus, test_sample_nexus)

    # checking if output directory exists, if it doesn't, creates the folder
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    # Load JSON for configuration
    # start with a fresh set of reduction parameters for every sample, because the reduction "pollutes"
    # the reduction parameters dictionary with additions not allowed by the schema
    reduction_input = json.loads(json_str)
    reduction_input["dataDirectories"] = nexus_dir
    reduction_input["configuration"]["outputDir"] = output_dir
    # reduction_input["sample"]["runNumber"] = source_sample_nexus
    reduction_input["sample"]["runNumber"] = test_sample_nexus
    reduction_input["sample"]["transmission"]["runNumber"] = samples_tran
    reduction_input["background"]["runNumber"] = backgrounds[0]
    reduction_input["background"]["transmission"]["runNumber"] = backgrounds_trans[0]
    reduction_input["outputFileName"] = sample_names[0]
    # always check after updating the parameters
    reduction_input = validate_reduction_parameters(reduction_input)
    loaded = load_all_files(reduction_input,
                            path=nexus_dir,
                            prefix=prefix)
    reduce_single_configuration(loaded, reduction_input)


def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix):

    unmatched_errors = ''

    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_test1_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        title_i = '{}: {}'.format(sample_name, title)
        try:
            compare_reduced_iq(output_log_file, gold_log_file, title_i, prefix)
        except AssertionError as unmatched_error:
            unmatched_errors = 'Testing output {} is different from gold result {}:\n{}' \
                               ''.format(output_log_file, gold_log_file, unmatched_error)
    # END-FOR

    # raise error for all
    if unmatched_errors != '':
        raise AssertionError(unmatched_errors)


def compare_reduced_iq(test_log_file, gold_log_file, title, prefix):
    """

    Parameters
    ----------
    test_log_file: str
        Absolute
    gold_log_file: str
    title: str
        plot title
    prefix: str
        file name prefix

    Returns
    -------

    """
    log_errors = list()

    for is_main_detector in [True, False]:
        assert os.path.exists(gold_log_file)
        vec_q_a, vec_i_a = get_iq1d(test_log_file, is_main=is_main_detector)
        vec_q_b, vec_i_b = get_iq1d(gold_log_file, is_main=is_main_detector)

        try:
            np.testing.assert_allclose(vec_q_a, vec_q_b)
            np.testing.assert_allclose(vec_i_a, vec_i_b)
            log_errors.append(None)
        except AssertionError as assert_err:
            log_errors.append(assert_err)
            from matplotlib import pyplot as plt
            if is_main_detector:
                flag = 'Main_detector'
            else:
                flag = 'Wing_detector'
            plt.cla()
            plt.plot(vec_q_a, vec_i_a, color='red', label='{} Corrected'.format(flag))
            plt.plot(vec_q_b, vec_i_b, color='black', label='{} Before being corrected'.format(flag))
            plt.yscale('log')
            plt.title(title)
            plt.legend()
            out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '_{}.png'.format(flag)
            plt.savefig(out_name)
    # END-FOR

    # Report
    if not (log_errors[0] is None and log_errors[1] is None):
        error_message = 'Main: {}; Wing: {}'.format(log_errors[0], log_errors[1])
        raise AssertionError(error_message)


def get_iq1d(log_file_name, is_main=True):
    """

    Parameters
    ----------
    log_file_name: str
        output log file's name
    is_main: bool
        for main or wing

    Returns
    -------

    """
    # Open file and entry
    assert os.path.exists(log_file_name)
    log_h5 = h5py.File(log_file_name, 'r')
    try:
        if is_main:
            iq1d_entry = log_h5['main_0']['I(Q)']
        else:
            iq1d_entry = log_h5['wing_0']['I(Q)']
    except KeyError:
        if is_main:
            iq1d_entry = log_h5['_slice_1']['main_0']['I(Q)']
        else:
            iq1d_entry = log_h5['_slice_1']['wing_0']['I(Q)']

    # Get data with a copy
    vec_q = np.copy(iq1d_entry['Q'].value)
    vec_i = np.copy(iq1d_entry['I'].value)

    # close file
    log_h5.close()

    return vec_q, vec_i


if __name__ == '__main__':
    pytest.main(__file__)
