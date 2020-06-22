import pytest
import numpy as np
import os
from drtsans.h5_buffer import parse_h5_entry
from drtsans.load import load_events
import h5py
from drtsans.mono.gpsans import (load_all_files, plot_reduction_output, reduce_single_configuration,
                                 reduction_parameters, update_reduction_parameters)
from drtsans.h5_buffer import FileNode, GroupNode, DataSetNode


def test_hello_world():
    assert 'Hello World!'


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


def next_test_copy_nexus(reference_dir, cleanfile):
    """Test creating event NeXus file, loading it and compare to the original event NeXus.

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
    nexus_h5.close()

    # Duplicate
    source_root.write(target_nexus)

    # Load source file to workspace
    source_ws = load_events(test_nexus_name, output_workspace='cg2_source')

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
    print('DEBUG Samples: {}'.format(samples[0]))
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


def copy_event_nexus_prototype01(source_nexus, target_nexus):

    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = parse_h5_entry(nexus_h5)
    # Duplicate
    source_root.write(target_nexus)
    # close HDF5
    nexus_h5.close()


def copy_event_nexus(source_nexus, target_nexus):
    """

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
                         # TODO '/entry/bank11_events',  # create node explicitly
                         '/entry/Software']

    # get children from entry node and duplicate except black list nodes
    for child_node in entry_node.children:
        if child_node.name not in level1_black_list:
            target_root_node.set_child(child_node)

    # set instrument node
    set_instrument_node(entry_node, target_entry_node)

    # set DAS logs
    set_das_log_node(entry_node, target_entry_node)

    # write
    target_root_node.write(target_nexus)

    # close original file
    nexus_h5.close()


def set_instrument_node(source_entry_node, target_entry_node):
    # now work with instrument
    source_instrument = source_entry_node.get_child('/entry/instrument')

    target_instrument = GroupNode(source_instrument.name)
    target_entry_node.set_child(target_instrument)
    target_instrument.add_attributes(source_instrument.attributes)
    # add all but not bank...
    for child in source_instrument.children:
        # skip all the bank*_events entries
        if child.name.split('/')[-1].startswith('bank'):
            print(f'skip node: {child.name}')
            continue   # not duplicating these entries
        target_instrument.set_child(child)


def set_das_log_node(source_entry_node, target_entry_node):
    target_logs_node = GroupNode('/entry/DASlogs')
    target_entry_node.set_child(target_logs_node)
    # add attribute
    target_logs_node.add_attributes({'NX_class': 'NXcollection'})

    # add sample logs
    source_logs_node = source_entry_node.get_child('/entry/DASlogs')

    # Specify white list
    logs_white_list = ['CG2:CS:SampleToSi',
                       'wavelength', 'wavelength_spread',
                       'sample_detector_distance',
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


def copy_event_nexus_prototype(source_nexus, target_nexus):
    """Try to use links and etc. to copy nexus

    Returns
    -------

    """

    # import source
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = parse_h5_entry(nexus_h5)

    # create a new file node
    target_root_node = FileNode()

    if True:
        # create an '/entry' node
        target_entry_node = GroupNode('/entry')
        target_root_node.set_child(target_entry_node)

        # set 'entry'
        entry_node = source_root.get_child('/entry')
        target_entry_node.add_attributes(entry_node.attributes)

        # define black_list
        black_list = ['/entry/user1',
                      '/entry/user2',
                      '/entry/user3',
                      '/entry/user4',
                      '/entry/user5',
                      '/entry/user6',
                      '/entry/user7',
                      '/entry/user8',
                      '/entry/instrument',
                      '/entry/Software',
                      '/entry/DASlogs']

        # get children from entry node
        for child_node in entry_node.children:
            if child_node.name not in black_list:
                target_root_node.set_child(child_node)

        # now work with instrument
        source_instrument = entry_node.get_child('/entry/instrument')

        target_instrument = GroupNode(source_instrument.name)
        target_entry_node.set_child(target_instrument)
        target_instrument.add_attributes(source_instrument.attributes)
        # add all but not bank...
        for child in source_instrument.children:
            if child.name.count('bank') > 0 and child.name.endswith('_events') > 0:
                print(child.name)
                continue   # not duplicating these entries
            target_instrument.set_child(child)

        # Das logs
        if True:
            target_logs_node = GroupNode('/entry/DASlogs')
            target_entry_node.set_child(target_logs_node)
            # add attribute
            target_logs_node.add_attributes({'NX_class': 'NXcollection'})

            # add sample logs
            source_logs_node = entry_node.get_child('/entry/DASlogs')

            # Test blacklist
            name_header_black_list = ['1K_Plate',
                                      'AllShutters',
                                      'Device',
                                      'Peltier',
                                      'ap',
                                      'guide',
                                      'CG2::SE:',
                                      'CG2::VS:',
                                      'CG2:Mot:', 'ILLF', 'Poly', 'coll', 'Cryo', 'He3', 'trap',
                                      'magr', 'p_',
                                      'Coll', 'CG2:SE', 'CG2:VS',
                                      'CG2:CS:',  # Partial black list
                                      'DLam', 'Lamb', 'Mag', 'Number', 'S', 'T', 'V',
                                      'a', 'b', 'c',
                                      'm', 'n', 'v', 'd',
                                      'f', 'p', 't',
                                      's', 'w'
                                      ]
            for child_log in source_logs_node.children:
                child_log_name = child_log.name.split('/')[-1]

                # Loop for black header
                on_black_list = False
                for black_header in name_header_black_list:
                    if child_log_name.startswith(black_header):
                        # skip the node with name's first several letters on black list
                        on_black_list = True
                        break

                # skip the node with name's first several letters on black list
                if on_black_list:
                    continue

                # add the node
                target_logs_node.set_child(child_log)

            # Test white list
            # logs_white_list = ['wavelength', 'wavelength_spread',
            #                    'CG2:CS:SampleToSi', 'sample_detector_distance',
            #                    'source_aperture_diameter', 'sample_aperture_diameter',
            #                    'proton_charge']
            logs_white_list = ['CG2:CS:SampleToSi',
                               # 'dcal', 'dcal_Readback',
                               # 'detector_trans',
                               'wavelength', 'wavelength_spread',
                               'sample_detector_distance',
                               'source_aperture_diameter', 'sample_aperture_diameter',
                               'detector_trans_Readback'] 
            for child_log in source_logs_node.children:
                child_log_name = child_log.name.split('/')[-1]

                if child_log_name in logs_white_list:
                    # only add nodes in white list
                    print(f'DEBUG SKIPPED: add DAS log {child_log.name}')
                    target_logs_node.set_child(child_log)
                else:
                    # target_logs_node.set_child(child_log)
                    continue

        else:
            # an alternative cannot-be-wrong case
            target_entry_node.set_child(entry_node.get_child('DASlogs'))

        # Print out the total number of sample logs
        print('[DEBUG] Number of DAS logs = {}'.format(len(target_logs_node.children)))

    else:
        # copy source entry
        source_entry = source_root.get_child('/entry')
        target_root_node.set_child(source_entry)

    # write
    target_root_node.write(target_nexus)

    # close original file
    nexus_h5.close()


def test_reduction(reference_dir):
    """Test reduction result between raw and generated event nexus file

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
    copy_event_nexus(source_nexus, target_nexus)

    # nexus_h5 = h5py.File(source_nexus, 'r')
    # source_root = parse_h5_entry(nexus_h5)
    # nexus_h5.close()
    # # Duplicate
    # source_root.write(target_nexus)

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
