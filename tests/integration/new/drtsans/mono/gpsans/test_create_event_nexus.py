"""
Integration test to create event nexus file
"""
import pytest
import numpy as np
import os
from drtsans.files.hdf5_rw import parse_h5_entry
from drtsans.load import load_events
import h5py
from drtsans.mono.gpsans import (load_all_files, plot_reduction_output, reduce_single_configuration,
                                 reduction_parameters, update_reduction_parameters)
from drtsans.files.hdf5_rw import FileNode, GroupNode
from drtsans.files.event_nexus_nodes import InstrumentNode, DasLogNode, BankNode
from drtsans.files.event_nexus_rw import convert_events_to_histogram, generate_events_from_histogram
from mantid.simpleapi import LoadEventNexus, SaveNexusProcessed


def test_step_by_step(reference_dir):
    # Generate a new event NeXus file (source)
    output_dir_src = '/tmp/nexus_src'
    if not os.path.exists(output_dir_src):
        os.mkdir(output_dir_src)

    # Generate a new event NeXus file
    output_dir = '/tmp/nexus'
    if not os.path.exists(output_dir):
        os.mkdir('/tmp/nexus')

    # Copy a nexus file
    test_nexus_name = 'CG2_9166.nxs.h5'
    source_nexus = os.path.join(reference_dir.new.gpsans, test_nexus_name)
    assert os.path.exists(source_nexus), f'Test data {source_nexus} does not exist'
    target_nexus = os.path.join(output_dir, 'CG2_9166.nxs.h5')

    # copy_event_nexus(source_nexus, target_nexus)
    generate_event_nexus(source_nexus, target_nexus)

    # Reduce
    loaded_raw_dict = reduce_data_step_by_step(source_nexus, reference_dir.new.gpsans, reference_dir.new.gpsans,
                                               output_dir_src, 'original')
    loaded_sample_ws = loaded_raw_dict.sample[0]
    print(f'sample ws: {loaded_sample_ws.name()}')
    print(f'sample ws type: {type(loaded_sample_ws)}')
    SaveNexusProcessed(InputWorkspace=loaded_sample_ws, Filename=os.path.join(output_dir_src, 'src_loaded.nxs'))
    print(f'spec nunber = {loaded_sample_ws.getNumberHistograms()}')

    # Reduce
    loaded_dup_dict = reduce_data_step_by_step(target_nexus, reference_dir.new.gpsans, reference_dir.new.gpsans,
                                               output_dir_src, 'duplicated')
    dup_sample_ws = loaded_dup_dict.sample[0]
    print(f'sample ws: {dup_sample_ws.name()}')
    SaveNexusProcessed(InputWorkspace=dup_sample_ws, Filename=os.path.join(output_dir, 'dup_loaded.nxs'))
    print(f'spec nunber = {loaded_sample_ws.getNumberHistograms()}')

    raw_x = loaded_sample_ws.extractX()
    dup_x = dup_sample_ws.extractX()
    diff_x = raw_x - dup_x
    print(np.where(np.abs(diff_x) > 0))
    print(diff_x[np.abs(diff_x) > 0]) 
    raw_y = loaded_sample_ws.extractY()
    dup_y = dup_sample_ws.extractY()
    diff_y = raw_y - dup_y
    print(np.where(np.abs(diff_y) > 0))
    print(diff_y[np.abs(diff_y) > 0]) 
    print(dup_y[np.abs(diff_y) > 0]) 
    np.testing.assert_allclose(raw_x, dup_x)
    np.testing.assert_allclose(raw_y, dup_y)

    raw_e = loaded_sample_ws.extractE()
    dup_e = dup_sample_ws.extractE()
    np.testing.assert_allclose(raw_e, dup_e)


def reduce_data_step_by_step(sample_nexus, data_dir, sens_dir, output_dir, prefix):

    sensitivity_file = os.path.join(sens_dir, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": 9177},
        "emptyTransmission": {"runNumber": 9177},
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

    # USER Input here with scan numbers etc.
    samples = [sample_nexus]  # ['9166']
    samples_trans = ['9178']
    sample_thick = ['0.1']
    bkgd = ['9165']
    bkgd_trans = ['9177']

    # Sample names for output
    sample_names = ["Al4"]

    # set output directory
    reduction_input["configuration"]["outputDir"] = output_dir
    # create output directory
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    specs = {
        "dataDirectories": data_dir,
        "sample": {"runNumber": samples[0],
                   "thickness": sample_thick[0],
                   "transmission": {"runNumber": samples_trans[0]}
                   },
        "background": {"runNumber": bkgd[0],
                       "transmission": {"runNumber": bkgd_trans[0]}
                       },
        "outputFileName": sample_names[0]
    }
    reduction_input = update_reduction_parameters(reduction_input, specs, validate=True)
    loaded = load_all_files(reduction_input, path=data_dir, prefix=prefix)

    print(f'LOADED: type = {type(loaded)}')

    # out = reduce_single_configuration(loaded, reduction_input)
    # plot_reduction_output(out, reduction_input, loglog=False)

    return loaded


def test_copy_h5_file(reference_dir, cleanfile):
    """Test duplicating an HDF5/NeXus.
    Verification is to load both original and duplicated file and compare to each other

    Test data: GPSANS run 9166

    Returns
    -------

    """
    # Get the source file
    test_nexus_name = 'CG2_9177.nxs.h5'
    source_nexus = os.path.join(reference_dir.new.gpsans, test_nexus_name)
    assert os.path.exists(source_nexus), f'Test data {source_nexus} does not exist'

    # Duplicate the source file to the temporary directory
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/nexus'
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir('/tmp/nexus')
    target_nexus = os.path.join(output_dir, 'CG2_9177.nxs.h5')

    # Load the source
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = parse_h5_entry(nexus_h5)

    # Duplicate
    source_root.write(target_nexus)
    nexus_h5.close()

    # Load source file to workspace
    source_ws = load_events(source_nexus, output_workspace='cg2_source')

    # Load the duplicated
    target_ws = load_events(target_nexus, output_workspace='cg2_duplicated')

    # Compare counts on each pixel
    source_y = source_ws.extractY()
    target_y = target_ws.extractY()
    np.testing.assert_allclose(source_y, target_y)

    # Compare pixels' positions
    num_hist = source_ws.getNumberHistograms()
    for iws in range(0, num_hist, 100):
        source_det_i_pos = source_ws.getInstrument().getDetector(iws).getPos()
        target_det_i_pos = target_ws.getInstrument().getDetector(iws).getPos()
        np.testing.assert_allclose(source_det_i_pos, target_det_i_pos,
                                   err_msg=f'Mismatch is detected at Detector {iws}')
    # Check source position
    source_moderator_pos = source_ws.getInstrument().getSource().getPos()
    target_moderator_pos = target_ws.getInstrument().getSource().getPos()
    np.testing.assert_allclose(source_moderator_pos, target_moderator_pos,
                               err_msg=f'Mismatch is detected at neutron source position')

    # Compare meta data
    assert len(source_ws.getRun().getProperties()) == len(target_ws.getRun().getProperties()), 'Meta data mismatch'


def reduce_gpsans_data(data_dir, reduction_input_common, output_dir, prefix, sample_nexus_path):
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
    # USER Input here with scan numbers etc.
    samples = [sample_nexus_path]  # ['9166']
    samples_trans = ['9178']
    sample_thick = ['0.1']
    bkgd = ['9165']
    bkgd_trans = ['9177']

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
    """Generate an event nexus file from source Nexus file

    This is the test case for various NeXus nodes

    Parameters
    ----------
    source_nexus
    target_nexus

    Returns
    -------

    """
    # import source
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = parse_h5_entry(nexus_h5)

    # create a new file node
    target_root_node = FileNode()

    # create an '/entry' node
    target_entry_node = GroupNode('/entry')
    target_root_node.set_child(target_entry_node)

    # set 'entry'
    entry_node = source_root.get_child('/entry')
    target_entry_node.add_attributes(entry_node.attributes)

    # define black_list under /entry directly
    level1_black_list = ['/entry/user1',
                         '/entry/user2',
                         '/entry/user3',
                         '/entry/user4',
                         '/entry/user5',
                         '/entry/user6',
                         '/entry/user7',
                         '/entry/user8',
                         '/entry/instrument',  # create node explicitly
                         '/entry/DASlogs',     # create node explicitly
                         '/entry/sample',
                         '/entry/entry_identifier',
                         '/entry/definition',
                         '/entry/total_uncounted_counts',
                         '/entry/bank9_events',  # create node explicitly
                         '/entry/Software']

    # get children from entry node and duplicate except black list nodes
    for child_node in entry_node.children:
        if child_node.name not in level1_black_list:
            target_root_node.set_child(child_node)

    # set instrument node
    set_instrument_node(nexus_h5, target_entry_node)

    # set DAS logs
    set_das_log_node(nexus_h5, entry_node, target_entry_node)

    # set Bank 9
    set_bank9_node(nexus_h5, target_entry_node)

    # write
    target_root_node.write(target_nexus)

    # close original file
    nexus_h5.close()


def set_bank9_node_exact_copy(source_h5, target_entry_node):
    """

    Parameters
    ----------
    source_h5: h5py._hl.files.File
        HDF5 file entry
    target_entry_node: GroupNode
        Target (output) group node for /entry/

    Returns
    -------

    """
    # Get a bank node
    bank9_entry = source_h5['/entry/bank9_events']
    event_ids = bank9_entry['event_id'][()]
    event_indexes = bank9_entry['event_index'][()]
    event_time_offsets = bank9_entry['event_time_offset'][()]
    event_time_zeros = bank9_entry['event_time_zero'][(())]
    run_start_time = bank9_entry['event_time_zero'].attrs['offset'].decode()

    # Create bank node for bank 9
    bank9_node = BankNode(name='/entry/bank9_events', bank_name='bank9')
    bank9_node.set_events(event_ids, event_indexes, event_time_offsets, run_start_time, event_time_zeros)

    # Link with parent
    target_entry_node.set_child(bank9_node)


def set_bank9_node(source_h5, target_entry_node):
    """Test writing bank 9 from histogram

    Parameters
    ----------
    source_h5: h5py._hl.files.File
        HDF5 file entry
    target_entry_node: GroupNode
        Target (output) group node for /entry/

    Returns
    -------

    """
    # Test with bank 9: retrieve information from bank 9
    bank9_entry = source_h5['/entry/bank9_events']
    bank9_histogram = convert_events_to_histogram(bank9_entry)
    run_start_time = bank9_entry['event_time_zero'].attrs['offset'].decode()

    # generate events
    nexus_events = generate_events_from_histogram(bank9_histogram, 10.)

    # Create bank node for bank 9
    bank9_node = BankNode(name='/entry/bank9_events', bank_name='bank9')
    bank9_node.set_events(nexus_events.event_id, nexus_events.event_index,
                          nexus_events.event_time_offset, run_start_time,
                          nexus_events.event_time_zero)

    # Link with parent
    target_entry_node.set_child(bank9_node)


def set_instrument_node(source_h5, target_entry_node):
    """Set instrument node

    Parameters
    ----------
    source_h5:  h5py._hl.files.File
        source entry node
    target_entry_node

    Returns
    -------

    """
    # IDF in XML
    xml_idf = source_h5['entry']['instrument']['instrument_xml']['data'][0]

    # Create new instrument node
    instrument_node = InstrumentNode()
    target_entry_node.set_child(instrument_node)

    # Set values
    instrument_node.set_idf(xml_idf, idf_type=b'text/xml', description=b'XML contents of the instrument IDF')
    instrument_node.set_instrument_info(target_station_number=1, beam_line=b'CG2', name=b'CG2', short_name=b'CG2')


def set_das_log_node(source_h5, source_entry_node, target_entry_node):
    """Set DAS log node in a mixed way

    Parameters
    ----------
    source_h5: H5 file
    source_entry_node: GroupNode
        source node
    target_entry_node: GroupNode
        target node

    Returns
    -------

    """
    target_logs_node = GroupNode('/entry/DASlogs')
    target_entry_node.set_child(target_logs_node)
    # add attribute
    target_logs_node.add_attributes({'NX_class': 'NXcollection'})

    # add sample logs
    source_logs_node = source_entry_node.get_child('/entry/DASlogs')

    # Specify white list
    logs_white_list = ['CG2:CS:SampleToSi',
                       'wavelength', 'wavelength_spread',
                       'source_aperture_diameter', 'sample_aperture_diameter',
                       'detector_trans_Readback']
    for child_log in source_logs_node.children:
        # remove HDF path from entry name
        child_log_name = child_log.name.split('/')[-1]

        if child_log_name in logs_white_list:
            # only add nodes in white list
            target_logs_node.set_child(child_log)
        else:
            # target_logs_node.set_child(child_log)
            continue

    # Add sample_detector_distance  manually
    set_sdd_node(target_logs_node, source_h5)


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
        print(error_message)
        raise AssertionError(error_message)


def test_reduction(reference_dir):
    """Test generate (partially copy) an event Nexus file by
    verifying reduction result between raw and generated event nexus file

    Testing is modified from mono.gpsans.test_overwrite_geometry_meta_data.test_no_overwrite()

    Returns
    -------

    """
    # Generate a new event NeXus file
    # TODO - in future it will be moved to a proper method in drtsans.generate_event_nexus
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/nexus'
    # cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir('/tmp/nexus')

    # Copy a nexus file
    test_nexus_name = 'CG2_9166.nxs.h5'
    source_nexus = os.path.join(reference_dir.new.gpsans, test_nexus_name)
    assert os.path.exists(source_nexus), f'Test data {source_nexus} does not exist'
    target_nexus = os.path.join(output_dir, 'CG2_9166.nxs.h5')

    # copy_event_nexus(source_nexus, target_nexus)
    generate_event_nexus(source_nexus, target_nexus)

    verify_histogram(source_nexus, target_nexus)

    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    # output_dir = mkdtemp(prefix='meta_overwrite_test1')
    # cleanfile(output_dir)
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": 9177},
        "emptyTransmission": {"runNumber": 9177},
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
                       sample_nexus_path=target_nexus)

    # Get result files
    sample_names = ["Al4"]
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test1/')

    # Verify results
    verify_reduction_results(sample_names, output_dir, gold_path,
                             title='Raw (No Overwriting)',  prefix='CG2MetaRaw')


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
    test_log_file
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


if __name__ == '__main__':
    pytest.main(__file__)
