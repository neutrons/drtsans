# Modified from gpsans/test_overwrite_geometry_meta_data.py
import pytest
import os
import json
import h5py
import numpy as np
from drtsans.mono.gpsans import load_all_files, reduce_single_configuration
from drtsans.mono.gpsans import plot_reduction_output
import time
from mantid.simpleapi import mtd


def reduce_gpsans_data(data_dir, json_file, output_dir, prefix, sample_detector_distance):
    """Standard reduction workflow

    Parameters
    ----------'
    data_dir
    json_file
    output_dir
    prefix: str
        prefix for all the workspaces loaded

    Returns
    -------

    """
    # USER Input here with scan numbers etc.
    samples = ['9176']
    samples_trans = ['9188']
    sample_thick = ['0.1'] * 3
    bkgd = ['9165']
    bkgd_trans = ['9177']

    # Sample names for output
    sample_names = ["PTMA-15"]

    # Import JSON
    with open(json_file) as f:
        reduction_input = json.load(f)

    # Set SDD
    reduction_input['configuration']['SampleDetectorDistance'] = '{}'.format(sample_detector_distance)

    # set output directory
    reduction_input["configuration"]["outputDir"] = output_dir
    # create output directory
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    start_time = time.time()
    for i in range(1):
        reduction_input["runNumber"] = samples[i]
        reduction_input["transmission"]["runNumber"] = samples_trans[i]
        reduction_input["background"]["runNumber"] = bkgd[i]
        reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[i]
        reduction_input["outputFilename"] = sample_names[i]
        reduction_input["thickness"] = sample_thick[i]
        loaded = load_all_files(reduction_input,
                                path=data_dir,
                                prefix=prefix)
        out = reduce_single_configuration(loaded, reduction_input)
        plot_reduction_output(out, reduction_input, loglog=False)

    end_time = time.time()
    print('Execution Time: {}'.format(end_time - start_time))


def test_multiple_sdd(reference_dir):
    """Test reduce 3 sets of data overwriting neither SampleToSi (distance) nor SampleDetectorDistance.

    This test case is provided by Lisa and verified by Lilin
    Location of original test: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/
    Test json:  /HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test1.json
    Verified result: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/test1/

    Returns
    -------

    """
    import shutil

    # Set test and run
    json_file = os.path.join('temp', 'gpsans_reduction_test1.json')
    assert os.path.exists(json_file), 'Test JSON {} cannot be accessed'.format(json_file)

    # Make final output by clearing out previously existing
    summary_path = 'Lilin_SDD_Test'
    if os.path.exists(summary_path):
        shutil.rmtree(summary_path)
    os.mkdir(summary_path)

    #  1m, 2m, 4m, 10m
    for tag, user_sdd in [('raw', ''),
                          ('1m', '1m')]:
        # Clean the data
        mtd.clear()

        # Set output
        output_dir = '/tmp/meta_overwrite_test_{}'.format(tag)
        reduce_gpsans_data(reference_dir.new.gpsans, json_file, output_dir, prefix='CG2Meta{}'.format(tag.upper()),
                           sample_detector_distance=user_sdd)

        # Copy files
        fig1d_name = os.path.join(output_dir, '1D/PTMA-15_1D.png')
        fig2d_name = os.path.join(output_dir, '2D/PTMA-15_2D.png')
        txt1d_name = os.path.join(output_dir, '1D/PTMA-15_1D.txt')
        for file_name in [fig1d_name, fig2d_name, txt1d_name]:
            name_parts = file_name.split('.')
            target_name = os.path.join(summary_path, name_parts[0] + '_' + tag + '_Before.' + name_parts[1])
            shutil.copy(file_name, target_name)
    # END-FOR


if __name__ == '__main__':
    pytest.main([__file__])
