"""
Integration test to create event nexus file
"""
import pytest
import numpy as np
import os
from drtsans.load import load_events
import h5py
from drtsans.mono.gpsans import (load_all_files, plot_reduction_output, reduce_single_configuration,
                                 reduction_parameters, update_reduction_parameters)
from drtsans.files.hdf5_rw import GroupNode, DataSetNode
from drtsans.files.event_nexus_nodes import InstrumentNode, DasLogNode, BankNode, MonitorNode
from drtsans.files.event_nexus_rw import generate_events_from_histogram
from drtsans.files.event_nexus_rw import generate_monitor_events_from_count
from drtsans.files.event_nexus_rw import init_event_nexus, parse_event_nexus, EventNeXusWriter
from drtsans.mono.convert_xml_to_nexus import EventNexusConverter
from mantid.simpleapi import LoadEventNexus, mtd, ConvertToMatrixWorkspace, LoadHFIRSANS
from tempfile import mkdtemp


def test_duplicate_event_nexus(reference_dir, cleanfile):
    """Test duplicating an HDF5/NeXus in 2 different approaches in order to verify EventNexusWriter

    Verification is to load both of the generated Event NeXus to do a comparison

    Test data: GPSANS run 9166

    Returns
    -------

    """
    # Get the source file
    source_nexus_file = 'CG2_9177.nxs.h5'
    source_nexus_file = os.path.join(reference_dir.new.gpsans, source_nexus_file)
    assert os.path.exists(source_nexus_file), f'Test data {source_nexus_file} does not exist'

    # Duplicate the source file to the temporary directory
    output_dir = mkdtemp(prefix='dupnexus')
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    prototype_dup_nexus = os.path.join(output_dir, 'CG2_9177_prototype.nxs.h5')
    product_dup_nexus = os.path.join(output_dir, 'CG2_9177_product.nxs.h5')

    # Duplicate with both approach
    generate_event_nexus_prototype(source_nexus_file, prototype_dup_nexus)
    generate_event_nexus(source_nexus_file, product_dup_nexus)

    # Load source file to workspace
    target_ws = load_events(product_dup_nexus, output_workspace='cg2_product', NumberOfBins=2)

    # Load the duplicated
    prototype_ws = load_events(prototype_dup_nexus, output_workspace='cg2_prototype', NumberOfBins=2)

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


def test_reduction(reference_dir, cleanfile):
    """Test generate (partially copy) an event Nexus file by
    verifying reduction result between raw and generated event nexus file

    Testing is modified from mono.gpsans.test_overwrite_geometry_meta_data.test_no_overwrite()

    Returns
    -------

    """
    # Generate a new event NeXus file
    # TODO - in future it will be moved to a proper method in drtsans.generate_event_nexus
    output_dir = mkdtemp(prefix='reducecg2nexus')
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Concert NeXus files
    nexus_file_dict = dict()
    for run_number in [9166, 9177, 9165, 9178]:
        # set source NeXus and target (testing) NeXus name
        test_nexus_name = f'CG2_{run_number}.nxs.h5'
        source_nexus = os.path.join(reference_dir.new.gpsans, test_nexus_name)
        assert os.path.exists(source_nexus), f'Test data {source_nexus} does not exist'
        target_nexus = os.path.join(output_dir, f'CG2_{run_number}.nxs.h5')
        # generate and verify
        generate_event_nexus(source_nexus, target_nexus)
        verify_histogram(source_nexus, target_nexus)
        # add to dictionary
        nexus_file_dict[run_number] = target_nexus

    # Set up reduction JSON
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": nexus_file_dict[9177]},
        "emptyTransmission": {"runNumber": nexus_file_dict[9177]},
        "configuration": {
            "outputDir": output_dir,
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-10,247-256'}"],
            "sensitivityFileName": sensitivity_file,
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": 40,
            "mmRadiusForTransmission": 40,
            "numQxQyBins": 180,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 180,
            "LogQBinsPerDecade": None,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'GPSANS', validate=False)  # add defaults and defer validation
    reduce_gpsans_data(reference_dir.new.gpsans, reduction_input, output_dir, prefix='CG2MetaRaw',
                       sample_nexus_path=nexus_file_dict[9166], sample_trans_path=nexus_file_dict[9178],
                       background_path=nexus_file_dict[9165], background_trans_path=nexus_file_dict[9177])

    # Get result files
    sample_names = ["Al4"]
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test1/')

    # Verify results
    verify_reduction_results(sample_names, output_dir, gold_path,
                             title='Raw (No Overwriting)',  prefix='CG2MetaRaw')


def reduce_gpsans_data(data_dir, reduction_input_common, output_dir, prefix, sample_nexus_path, sample_trans_path,
                       background_path, background_trans_path):
    """Standard reduction workflow

    Parameters
    ----------'
    data_dir
    reduction_input_common: dict
        reduction parameters common to all samples
    output_dir
    prefix: str
        prefix for all the workspaces loaded

    Returns
    -------

    """
    # sample_trans_file = None):
    # USER Input here with scan numbers etc.
    samples = [sample_nexus_path]  # ['9166']
    samples_trans = [sample_trans_path]  # ['9178']
    sample_thick = ['0.1']
    bkgd = [background_path]  # 9165
    bkgd_trans = [background_trans_path]  # ['9177']

    # Sample names for output
    sample_names = ["Al4"]

    # set output directory
    reduction_input_common["configuration"]["outputDir"] = output_dir
    # create output directory
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    for i in range(len(samples)):
        specs = {
            "dataDirectories": data_dir,
            "sample": {"runNumber": samples[i],
                       "thickness": sample_thick[i],
                       "transmission": {"runNumber": samples_trans[i]}
                       },
            "background": {"runNumber": bkgd[i],
                           "transmission": {"runNumber": bkgd_trans[i]}
                           },
            "outputFileName": sample_names[i]
        }
        reduction_input = update_reduction_parameters(reduction_input_common, specs, validate=True)
        loaded = load_all_files(reduction_input, path=data_dir, prefix=prefix)
        out = reduce_single_configuration(loaded, reduction_input)
        plot_reduction_output(out, reduction_input, loglog=False)


def generate_event_nexus(source_nexus, target_nexus):
    """Generate event NeXus properly

    Parameters
    ----------
    source_nexus
    target_nexus

    Returns
    -------

    """
    # Import essential experimental data from source event nexus file
    nexus_contents = parse_event_nexus(source_nexus, num_banks=48)
    # Generate event nexus writer
    event_nexus_writer = EventNeXusWriter(beam_line='CG2', instrument_name='CG2')

    # set instrument
    event_nexus_writer.set_instrument_info(48,  nexus_contents[0])

    # set counts
    for bank_id in range(1, 48 + 1):
        event_nexus_writer.set_bank_histogram(bank_id, nexus_contents[1][bank_id])

    # set meta
    for das_log in nexus_contents[5].values():
        event_nexus_writer.set_meta_data(das_log)

    # time
    start_time = nexus_contents[3]
    end_time = nexus_contents[4]

    # Write file
    event_nexus_writer.generate_event_nexus(target_nexus, start_time, end_time, nexus_contents[2])


def generate_event_nexus_prototype(source_nexus, target_nexus):
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
    nexus_contents = parse_event_nexus(source_nexus, num_banks=48)

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
        ('/entry/end_time', nexus_contents[4]),
        ('/entry/duration', 'only used by biosans')
    ]
    for child_node_name, child_value in entry_level_white_list:
        child_node = DataSetNode(child_node_name)
        child_node.set_string_value(child_value)
        target_entry_node.set_child(child_node)

    # set Bank 1 - 48
    max_pulse_time_array = None
    for bank_id in range(1, 48 + 1):
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
        bank ID (from 1 to 48)
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
        run_start_time = run_start_time.decode()
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
    instrument_node.set_instrument_info(target_station_number=1, beam_line=b'CG2', name=b'CG2', short_name=b'CG2')


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
                                  target=b'/entry/DASlogs/CG2:CS:SampleToDetRBV')

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
    if error_message != '':
        raise AssertionError(error_message)


def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix):

    unmatched_errors = ''

    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
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
    """Compare I(Q) from reduced file and gold file

    Parameters
    ----------
    test_log_file: str
        name of drtsans result log file
    gold_log_file
    title: str
        title of output figure
    prefix: str
        prefix of output file

    Returns
    -------

    """
    # Plot main
    test_q_vec, test_intensity_vec = get_iq1d(test_log_file)
    gold_q_vec, gold_intensity_vec = get_iq1d(gold_log_file)

    # Verify result
    try:
        np.testing.assert_allclose(test_q_vec, test_q_vec, atol=1E-4)
        np.testing.assert_allclose(test_intensity_vec, gold_intensity_vec, atol=1E-7)
    except AssertionError as assert_err:
        from matplotlib import pyplot as plt
        plt.cla()
        plt.plot(test_q_vec, test_intensity_vec, color='red', label='Corrected')
        plt.plot(gold_q_vec, gold_intensity_vec, color='black', label='Before being corrected')
        plt.legend()
        plt.title(title)
        plt.yscale('log')
        out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '.png'
        plt.savefig(out_name)

        raise assert_err


def get_iq1d(log_file_name):
    """Get I(Q) from output SANS log file

    Parameters
    ----------
    log_file_name: str
        log file name

    Returns
    -------
    tuple
        numpy 1D array for Q, numpy 1D array for intensity

    """
    # Open file and entry
    log_h5 = h5py.File(log_file_name, 'r')

    if '_slice_1' in log_h5:
        data_entry = log_h5['_slice_1']['main']
    else:
        data_entry = log_h5['main']

    # Get data
    iq1d_entry = data_entry['I(Q)']

    # Get data with a copy
    vec_q = np.copy(iq1d_entry['Q'][()])
    vec_i = np.copy(iq1d_entry['I'][()])

    # close file
    log_h5.close()

    return vec_q, vec_i


def test_convert_spice_to_nexus(reference_dir, cleanfile):
    """Test to convert SPICE to NeXus

    Parameters
    ----------
    reference_dir
    cleanfile

    Returns
    -------

    """
    # Specify the test data
    spice_data_file = os.path.join(reference_dir.new.gpsans, 'CG2_exp315_scan0005_0060.xml')
    template_nexus_file = os.path.join(reference_dir.new.gpsans, 'CG2_9177.nxs.h5')
    assert os.path.exists(spice_data_file)
    assert os.path.exists(template_nexus_file)

    output_dir = mkdtemp(prefix='spice2nexus')
    cleanfile(output_dir)

    # Convert from SPICE to event Nexus
    out_nexus_file = os.path.join(output_dir, 'CG2_31500050060.nxs.h5')

    # init convert
    converter = EventNexusConverter('CG2', 'CG2')
    converter.load_idf(template_nexus_file)
    converter.load_sans_xml(spice_data_file)
    converter.generate_event_nexus(out_nexus_file, num_banks=48)

    # Check
    os.path.exists(out_nexus_file)

    # Check instrument node against the original one
    test_nexus_h5 = h5py.File(out_nexus_file, 'r')
    test_idf = test_nexus_h5['entry']['instrument']['instrument_xml']['data'][0]
    expected_nexus_h5 = h5py.File(template_nexus_file, 'r')
    expected_idf = expected_nexus_h5['entry']['instrument']['instrument_xml']['data'][0]
    assert test_idf == expected_idf
    test_nexus_h5.close()
    expected_nexus_h5.close()

    # Load
    test_ws_name = 'TestSpice2Nexus315560'
    LoadEventNexus(Filename=out_nexus_file, OutputWorkspace=test_ws_name,
                   NumberOfBins=1, LoadNexusInstrumentXML=True)
    ConvertToMatrixWorkspace(InputWorkspace=test_ws_name, OutputWorkspace=test_ws_name)
    test_nexus_ws = mtd[test_ws_name]

    # Load template event nexus
    LoadEventNexus(Filename=template_nexus_file, OutputWorkspace='cg3template',
                   NumberOfBins=1, LoadNexusInstrumentXML=True)
    template_ws = mtd['cg3template']

    # Check number of histograms
    assert test_nexus_ws.getNumberHistograms() == template_ws.getNumberHistograms()

    # Compare units of required DAS logs
    for das_log_name in ['CG2:CS:SampleToSi', 'wavelength', 'wavelength_spread', 'source_aperture_diameter',
                         'sample_aperture_diameter', 'detector_trans_Readback', 'sample_detector_distance',
                         'detector_trans_Readback']:
        template_unit = template_ws.run().getProperty(das_log_name).units
        test_unit = test_nexus_ws.run().getProperty(das_log_name).units
        assert template_unit == test_unit, f'DAS log {das_log_name} unit does not match'

    # Check instrument by comparing pixel position
    # Run 9711: detector_trans_Readback = 0.002 mm (to negative X direction)
    # Exp315 Scan 5  Run 60: detector trans = 0.001 mm
    # Thus all pixels of from-SPICE data shall have a postive 1 mm shift
    # Both data have different SDD.  Thus all the pixels will have a constant shift along Z direction
    diff_x = 0.001
    diff_z_list = list()
    for iws in range(0, template_ws.getNumberHistograms(), 10):
        test_pixel_pos = test_nexus_ws.getDetector(iws).getPos()
        expected_pixel_pos = template_ws.getDetector(iws).getPos()
        # constant difference at x
        assert test_pixel_pos[0] - diff_x == pytest.approx(expected_pixel_pos[0], abs=1e-7)
        # y shall be exactly same
        assert test_pixel_pos[1] == pytest.approx(expected_pixel_pos[1], abs=1e-7)
        # z shall have constant difference
        diff_z_list.append(test_pixel_pos[2] - expected_pixel_pos[2])
    # shift along Z-axis shall be a constant
    assert np.array(diff_z_list).std() < 1E-12

    # Load original SPICE file
    spice_ws_name = os.path.basename(spice_data_file).split('.')[0]
    spice_ws_name = f'CG3IntTestSpice_{spice_ws_name}'
    LoadHFIRSANS(Filename=spice_data_file, OutputWorkspace=spice_ws_name)
    spice_ws = mtd[spice_ws_name]

    # compare histograms
    for iws in range(0, test_nexus_ws.getNumberHistograms()):
        assert test_nexus_ws.readY(iws)[0] == pytest.approx(spice_ws.readY(iws + 2)[0], abs=1E-3)

    # compare DAS logs (partial)
    for log_name in ['wavelength', 'source_aperture_diameter', 'sample_aperture_diameter']:
        nexus_log_value = test_nexus_ws.run().getProperty(log_name).value.mean()
        spice_log_value = spice_ws.run().getProperty(log_name).value
        assert nexus_log_value == pytest.approx(spice_log_value, 1e-7)


if __name__ == '__main__':
    pytest.main(__file__)
