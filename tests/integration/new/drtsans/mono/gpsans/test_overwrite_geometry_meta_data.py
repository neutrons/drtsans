# Integration test for overwriting instrument geometry related meta data for GP-SANS
# Tests are
import pytest
import os
import json
import h5py
import numpy as np
from drtsans.mono.gpsans import load_all_files, reduce_single_configuration
import time
from mantid.simpleapi import mtd


def reduce_gpsans_data(json_file, output_dir):
    """Standard reduction workflow

    Parameters
    ----------
    json_file
    output_dir

    Returns
    -------

    """
    # Clear the existing workspaces to force reloading data for various geometry setup
    clean_workspaces()

    # USER Input here with scan numbers etc.
    samples = ['9166', '9167', '9176']
    samples_trans = ['9178', '9179', '9188']
    sample_thick = ['0.1'] * 3
    bkgd = ['9165', '9165', '9165']
    bkgd_trans = ['9177', '9177', '9177']

    # Sample names for output
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]

    # Import JSON
    with open(json_file) as f:
        reduction_input = json.load(f)

    # set output directory
    reduction_input["configuration"]["outputDir"] = output_dir
    # create output directory
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    start_time = time.time()
    for i in range(len(samples)):
        reduction_input["runNumber"] = samples[i]
        reduction_input["transmission"]["runNumber"] = samples_trans[i]
        reduction_input["background"]["runNumber"] = bkgd[i]
        reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[i]
        reduction_input["outputFilename"] = sample_names[i]
        reduction_input["thickness"] = sample_thick[i]
        loaded = load_all_files(reduction_input,
                                path='/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/gpsans/data')
        out = reduce_single_configuration(loaded, reduction_input)
        assert out
        # plot_reduction_output(out, reduction_input, loglog=False)

    end_time = time.time()
    print('Execution Time: {}'.format(end_time - start_time))


def clean_workspaces():
    """Clean all the workspaces in AnalysisDataService

    Returns
    -------

    """
    # workspace_names = mtd.getObjectNames()
    # for ws_name in workspace_names:
    #     mtd.remove(ws_name)
    mtd.clear()


def get_iq1d(log_file_name):
    """

    Parameters
    ----------
    log_file_name

    Returns
    -------

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
    vec_q = np.copy(iq1d_entry['Q'].value)
    vec_i = np.copy(iq1d_entry['I'].value)

    # close file
    log_h5.close()

    return vec_q, vec_i


def compare_reduced_iq(test_log_file, gold_log_file):
    """Compare I(Q) from reduced file and gold file

    Parameters
    ----------
    test_log_file
    gold_log_file

    Returns
    -------

    """
    # Plot main
    test_q_vec, test_intensity_vec = get_iq1d(test_log_file)
    gold_q_vec, gold_intensity_vec = get_iq1d(gold_log_file)

    # Verify result
    np.testing.assert_allclose(test_q_vec, test_q_vec, atol=1E-4)
    np.testing.assert_allclose(test_intensity_vec, gold_intensity_vec, atol=1E-7)


def verify_reduction_results(sample_names, output_dir, gold_path):
    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        compare_reduced_iq(output_log_file, gold_log_file)


def test_no_overwrite():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set test and run
    json_file =\
        '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/gpsans_reduction_test1.json'
    output_dir = '/tmp/meta_overwrite_test1'
    reduce_gpsans_data(json_file, output_dir)

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/test1/'

    # Verify results
    verify_reduction_results(sample_names, output_dir, gold_path)

    # for sample_name in sample_names:
    # # output log file name
    # output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
    # assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
    # # gold file
    # gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
    # assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
    # # compare
    # compare_reduced_iq(output_log_file, gold_log_file)


def test_overwrite_sample2si():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set test and run: sample to silicon window is changed 94 mm
    json_file = \
        '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/gpsans_reduction_test2.json'
    output_dir = '/tmp/meta_overwrite_test2'
    reduce_gpsans_data(json_file, output_dir)

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]

    # Verify results
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/test2/'
    verify_reduction_results(sample_names, output_dir, gold_path)
    # for sample_name in sample_names:
    #     # output log file name
    #     output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
    #     assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
    #     # gold file
    #     gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
    #     assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
    #     # compare
    #     compare_reduced_iq(output_log_file, gold_log_file)


def skip_test_overwrite_sdd():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set test and run: sample to detector distance is changed to 40 meter
    json_file = \
        '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/gpsans_reduction_test3.json'
    output_dir = '/tmp/meta_overwrite_test3'
    reduce_gpsans_data(json_file, output_dir)

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    output_log_files = [os.path.join(output_dir, '{}_reduction_log.hdf'.format(sn)) for sn in sample_names]
    for output_file_path in output_log_files:
        assert os.path.exists(output_file_path), 'Output {} cannot be found'.format(output_file_path)

    # Verify results
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/test3/'
    verify_reduction_results(sample_names, output_dir, gold_path)


def skip_test_overwrite_both():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set test and run: sample to silicon window to 94 mm and sample to detector distance to 15 meter
    json_file = \
        '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/gpsans_reduction_test4.json'
    output_dir = '/tmp/meta_overwrite_test4'
    reduce_gpsans_data(json_file, output_dir)

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    output_log_files = [os.path.join(output_dir, '{}_reduction_log.hdf'.format(sn)) for sn in sample_names]
    for output_file_path in output_log_files:
        assert os.path.exists(output_file_path), 'Output {} cannot be found'.format(output_file_path)

    # Verify results
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/test4/'
    verify_reduction_results(sample_names, output_dir, gold_path)


if __name__ == '__main__':
    pytest.main([__file__])
