# Integration test for overwriting instrument geometry related meta data for GP-SANS
import pytest
from tempfile import mkdtemp
import time
import os
from drtsans.files.log_h5_reader import verify_cg2_reduction_results

from drtsans.mono.gpsans import (load_all_files, plot_reduction_output, reduce_single_configuration,
                                 reduction_parameters, update_reduction_parameters)


def reduce_gpsans_data(data_dir, reduction_input_common, output_dir, prefix):
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
    samples = ['9166', '9167', '9176']
    samples_trans = ['9178', '9179', '9188']
    sample_thick = ['0.1'] * 3
    bkgd = ['9165', '9165', '9165']
    bkgd_trans = ['9177', '9177', '9177']

    # Sample names for output
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]

    # set output directory
    reduction_input_common["configuration"]["outputDir"] = output_dir
    # create output directory
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    start_time = time.time()
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

    end_time = time.time()
    print('Execution Time: {}'.format(end_time - start_time))


# def get_iq1d(log_file_name):
#     """Get I(Q) from output SANS log file
#
#     Parameters
#     ----------
#     log_file_name: str
#         log file name
#
#     Returns
#     -------
#     tuple
#         numpy 1D array for Q, numpy 1D array for intensity
#
#     """
#     # Open file and entry
#     log_h5 = h5py.File(log_file_name, 'r')
#
#     if '_slice_1' in log_h5:
#         data_entry = log_h5['_slice_1']['main']
#     else:
#         data_entry = log_h5['main']
#
#     # Get data
#     iq1d_entry = data_entry['I(Q)']
#
#     # Get data with a copy
#     vec_q = np.copy(iq1d_entry['Q'].value)
#     vec_i = np.copy(iq1d_entry['I'].value)
#
#     # close file
#     log_h5.close()
#
#     return vec_q, vec_i
#
#
# def compare_reduced_iq(test_log_file, gold_log_file, title, prefix):
#     """Compare I(Q) from reduced file and gold file
#
#     Parameters
#     ----------
#     test_log_file
#     gold_log_file
#     title: str
#         title of output figure
#     prefix: str
#         prefix of output file
#
#     Returns
#     -------
#
#     """
#     # Plot main
#     test_q_vec, test_intensity_vec = get_iq1d(test_log_file)
#     gold_q_vec, gold_intensity_vec = get_iq1d(gold_log_file)
#
#     # Verify result
#     try:
#         np.testing.assert_allclose(test_q_vec, test_q_vec, atol=1E-4)
#         np.testing.assert_allclose(test_intensity_vec, gold_intensity_vec, atol=1E-7)
#     except AssertionError as assert_err:
#         from matplotlib import pyplot as plt
#         plt.cla()
#         plt.plot(test_q_vec, test_intensity_vec, color='red', label='Corrected')
#         plt.plot(gold_q_vec, gold_intensity_vec, color='black', label='Before being corrected')
#         plt.legend()
#         plt.title(title)
#         plt.yscale('log')
#         out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '.png'
#         plt.savefig(out_name)
#
#         raise assert_err
#
#
# def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix):
#
#     unmatched_errors = ''
#
#     for sample_name in sample_names:
#         # output log file name
#         output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
#         assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
#         # gold file
#         gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
#         assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
#         # compare
#         title_i = '{}: {}'.format(sample_name, title)
#         try:
#             compare_reduced_iq(output_log_file, gold_log_file, title_i, prefix)
#         except AssertionError as unmatched_error:
#             unmatched_errors = 'Testing output {} is different from gold result {}:\n{}' \
#                                ''.format(output_log_file, gold_log_file, unmatched_error)
#     # END-FOR
#
#     # raise error for all
#     if unmatched_errors != '':
#         raise AssertionError(unmatched_errors)


# dev - Wenduo Zhou <wzz@ornl.gov>
# SME - Debeer-Schmitt, Lisa M. debeerschmlm@ornl.gov, He, Lilin <hel3@ornl.gov>
def test_no_overwrite(reference_dir, cleanfile):
    """Test reduce 3 sets of data overwriting neither SampleToSi (distance) nor SampleDetectorDistance.

    This test case is provided by Lisa and verified by Lilin
    Location of original test: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/
    Test json:  /HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test1.json
    Verified result: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/test1/

    Returns
    -------

    """
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    output_dir = mkdtemp(prefix='meta_overwrite_test1')
    cleanfile(output_dir)
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
    reduce_gpsans_data(reference_dir.new.gpsans, reduction_input, output_dir, prefix='CG2MetaRaw')

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test1/')

    # Verify results
    verify_cg2_reduction_results(sample_names, output_dir, gold_path,
                                 title='Raw (No Overwriting)', prefix='CG2MetaRaw')


# dev - Wenduo Zhou <wzz@ornl.gov>
# SME - Debeer-Schmitt, Lisa M. debeerschmlm@ornl.gov, He, Lilin <hel3@ornl.gov>
def test_overwrite_sample2si(reference_dir, cleanfile):
    """Test reduce 3 sets of data overwriting SampleToSi (distance) but not SampleDetectorDistance.
    Sample to detector distance will be changed accordingly.

    - Overwrite SampleToSi (distance) to 94 mm.

    This test case is provided by Lisa and verified by Lilin
    Location of original test: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/
    Test json:  /HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test2.json
    Verified result: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/test2/

    Returns
    -------

    """
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    output_dir = mkdtemp(prefix='meta_overwrite_test2')
    cleanfile(output_dir)
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": 9177},
        "emptyTransmission": {"runNumber": 9177},
        "configuration": {
            "outputDir": output_dir,
            "sampleToSi": 94,
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-10,247-256'}"],
            "sensitivityFileName": sensitivity_file,
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": 40,
            "mmRadiusForTransmission": 40,
            "numQxQyBins": 150,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 150,
            "LogQBinsPerDecade": None,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'GPSANS', validate=False)  # add defaults and defer validation
    reduce_gpsans_data(reference_dir.new.gpsans, reduction_input, output_dir, 'CG2MetaSWD')

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]

    # Verify results
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test2/')
    verify_cg2_reduction_results(sample_names, output_dir, gold_path,
                                 title='Overwrite SampleToSi to 94mm', prefix='CG2MetaSWD')


# dev - Wenduo Zhou <wzz@ornl.gov>
# SME - Debeer-Schmitt, Lisa M. debeerschmlm@ornl.gov, He, Lilin <hel3@ornl.gov>
def test_overwrite_sdd(reference_dir, cleanfile):
    """Test reduce 3 sets of data overwriting SampleDetectorDistance but not SampleDetectorDistance

    - Overwrite DetectorToSample (distance) to 40 meter

    This test case is provided by Lisa and verified by Lilin
    Location of original test: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/
    Test json:  /HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test3.json
    Verified result: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/test3/

    Returns
    -------

    """
    # Set test and run: sample to detector distance is changed to 40 meter
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    output_dir = mkdtemp(prefix='meta_overwrite_test3')
    cleanfile(output_dir)
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": 9177},
        "emptyTransmission": {"runNumber": 9177},
        "configuration": {
            "outputDir": output_dir,
            "sampleDetectorDistance": 40,
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-10,247-256'}"],
            "sensitivityFileName": sensitivity_file,
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": 40,
            "mmRadiusForTransmission": 40,
            "numQxQyBins": 150,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 150,
            "LogQBinsPerDecade": None,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'GPSANS', validate=False)  # add defaults and defer validation
    reduce_gpsans_data(reference_dir.new.gpsans, reduction_input, output_dir, 'CG2MetaSDD')

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    output_log_files = [os.path.join(output_dir, '{}_reduction_log.hdf'.format(sn)) for sn in sample_names]
    for output_file_path in output_log_files:
        assert os.path.exists(output_file_path), 'Output {} cannot be found'.format(output_file_path)

    # Verify results
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test3/')
    verify_cg2_reduction_results(sample_names, output_dir, gold_path,
                                 title='Overwrite DetectorSampleDistance to 40 meter',
                                 prefix='CG2MetaSDD')


# dev - Wenduo Zhou <wzz@ornl.gov>
# SME - Debeer-Schmitt, Lisa M. debeerschmlm@ornl.gov, He, Lilin <hel3@ornl.gov>
def test_overwrite_both(reference_dir, cleanfile):
    """Test reduce 3 sets of data overwriting both SampleToSi (distance) and SampleDetectorDistance

    - Overwrite SampleToSi (distance) to 200 mm.
    - Overwrite DetectorToSample (distance) to 30 meter

    This test case is provided by Lisa and verified by Lilin
    Location of original test: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/
    Test json:  /HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test4.json
    Verified result: /HFIR/CG2/shared/UserAcceptance/overwrite_meta/test4/

    Returns
    -------

    """
    # Set test and run: sample to silicon window to 94 mm and sample to detector distance to 15 meter
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/sens_c486_noBar.nxs')
    output_dir = mkdtemp(prefix='meta_overwrite_test4')
    cleanfile(output_dir)
    specs = {
        "iptsNumber": 21981,
        "beamCenter": {"runNumber": 9177},
        "emptyTransmission": {"runNumber": 9177},
        "configuration": {
            "outputDir": output_dir,
            "sampleToSi": 200,
            "sampleDetectorDistance": 30,
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-10,247-256'}"],
            "sensitivityFileName": sensitivity_file,
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": 40,
            "mmRadiusForTransmission": 40,
            "numQxQyBins": 150,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 150,
            "LogQBinsPerDecade": None,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'GPSANS', validate=False)  # add defaults and defer validation
    reduce_gpsans_data(reference_dir.new.gpsans, reduction_input, output_dir, 'CG2MetaBoth')

    # Get result files
    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    output_log_files = [os.path.join(output_dir, '{}_reduction_log.hdf'.format(sn)) for sn in sample_names]
    for output_file_path in output_log_files:
        assert os.path.exists(output_file_path), 'Output {} cannot be found'.format(output_file_path)

    # Verify results
    gold_path = os.path.join(reference_dir.new.gpsans, 'overwrite_gold_04282020/test4/')
    verify_cg2_reduction_results(sample_names, output_dir, gold_path,
                                 title='Overwrite DetectorSampleDistance to 30 meter, SampleToSi to 200 mm',
                                 prefix='CG2MetaBoth')


if __name__ == '__main__':
    pytest.main([__file__])
