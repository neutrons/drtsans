# Integration test for overwriting instrument geometry related meta data for BIO-SANS
# From round 3 4822
# Tests are
import pytest
import numpy as np
import json
import os
import h5py
from drtsans.mono.biosans import load_all_files, reduce_single_configuration
import time
from mantid.simpleapi import mtd


def test_no_overwrite():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set up test
    json_file = generate_testing_json(None, None)
    output_dir = '/tmp/meta_overwrite_bio_test1/'

    # Run
    reduce_biosans_data(json_file, output_dir)

    # Get result files
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/test1/'
    # gold_path = '/HFIR/CG3/shared/UserAcceptance/override_round3'

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, 'test1')


def test_overwrite_both_minor():
    """Test reduce 3 sets of data with minor overwriting

    Returns
    -------

    """
    # Set up test
    json_file = generate_testing_json(61, 6.9)
    output_dir = '/tmp/meta_overwrite_bio_test1a/'

    # Run
    reduce_biosans_data(json_file, output_dir)

    # Get result files
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/test1a/'

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, 'test1a')


def test_overwrite_both_major():
    """Test reduce 3 sets of data with minor overwriting

    Returns
    -------

    """
    # Set up test
    json_file = generate_testing_json(200, 14)
    output_dir = '/tmp/meta_overwrite_bio_test4/'

    # Run
    reduce_biosans_data(json_file, output_dir)

    # Get result files
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/test4/'

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, 'test4')


def test_overwrite_sample_to_si():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set up test
    json_file = generate_testing_json(7000, None)
    output_dir = '/tmp/meta_overwrite_bio_test2/'

    # Run
    reduce_biosans_data(json_file, output_dir)

    # Get result files
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/test2/'

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, 'test2')


def test_overwrite_sample_to_detector():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set up test
    json_file = generate_testing_json(None, 14)
    output_dir = '/tmp/meta_overwrite_bio_test3/'

    # Run
    reduce_biosans_data(json_file, output_dir)

    # Get result files
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    gold_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/test3/'

    # Verify
    verify_reduction_results(sample_names, output_dir, gold_path, 'test3')


def reduce_biosans_data(json_str, output_dir):
    """

    Parameters
    ----------
    json_str
    output_dir

    Returns
    -------

    """
    # Clear workspaces in memory
    mtd.clear()

    # Set up (testing) runs
    sample_names = ['csmb_ecoli1h_n2', 'insect1hTime_n2']
    samples = ['5709', '5712']
    samples_trans = samples
    backgrounds = ['5715', '5715']
    backgrounds_trans = backgrounds

    # Test data path
    nexus_path = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/biosans/data'

    # Load JSON for configuration
    reduction_input = json.loads(json_str)

    # chekcing if output directory exists, if it doesn't, creates the folder
    reduction_input["configuration"]["outputDir"] = output_dir
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    start_time = time.time()
    for i in range(len(samples)):
        reduction_input["runNumber"] = samples[i]
        reduction_input["transmission"]["runNumber"] = samples_trans[i]
        reduction_input["background"]["runNumber"] = backgrounds[i]
        reduction_input["background"]["transmission"]["runNumber"] = backgrounds_trans[i]
        reduction_input["outputFilename"] = sample_names[i]
        loaded = load_all_files(reduction_input,
                                path=nexus_path)
        out = reduce_single_configuration(loaded, reduction_input)
        assert out

    end_time = time.time()
    print(end_time - start_time)


def generate_testing_json(sample_to_si_window_distance, sample_to_detector_distance):
    """Generating testing JSON

    Parameters
    ----------
    sample_to_si_window_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit millimeter
    sample_to_detector_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit meter

    Returns
    -------
    str
        JSON string

    """
    json_str = """ {
    "instrumentName": "CG3",
    "iptsNumber": "24740",
    "runNumber": "4822",
    "thickness": "0.1",
    "outputFilename": "CG3_4822",
    "transmission": {
        "runNumber": "4822",
        "value": ""
    },
    "background": {
        "runNumber": "4821",
        "transmission": {
            "runNumber": "4821",
            "value": ""
        }
    },
    "beamCenter": {
        "runNumber": "1322"
    },
    "emptyTrans": {
        "runNumber": "5705"
    },
    "configuration": {
        "outputDir": "/HFIR/CG3/shared/UserAcceptance/override/test1",
        "sampleApertureSize": "14",
        "sourceApertureDiameter": "",
        "maskFileName": "",
        "useMaskFileName": false,
        "useDefaultMask": true,
        "DefaultMask":["{'Pixel':'1-12,244-256'}", "{'Bank':'21-24,45-48'}"],
        "useBlockedBeam": false,
        "BlockBeamFileName":"",
        "useDarkFileName": true,
        "darkMainFileName": "CG3_1383.nxs",
        "darkWingFileName": "CG3_1383.nxs",
        "useSensitivityFileName": true,
        "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle486/sens_f4829m7p0_TDC_SAC.h5",
        "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle486/sens_f4835w3p2_TDC_SAC.h5",
        "UseBarScan": false,
        "BarScanMainFileName":"",
        "BarScanWingFileName":"",
        "absoluteScaleMethod":"standard",
        "DBScalingBeamRadius": "",
        "StandardAbsoluteScale": "0.0055e-8",
        "normalization": "Monitor",
        "sampleOffset": "",
        "useSampleOffset": false,
        "useDetectorTubeType": true,
        "useSolidAngleCorrection": true,
        "useThetaDepTransCorrection": true,
        "mmRadiusForTransmission": "",
        "numMainQxQyBins": "100",
        "numWingQxQyBins": "100",
        "1DQbinType": "scalar",
        "QbinType": "log",
        "LogQBinsEvenDecade": false,
        "LogQBinsPerDecadeMain":20,
        "LogQBinsPerDecadeWing": 25,
        "WedgeMinAngles": "-30, 60",
        "WedgeMaxAngles": "30, 120",
        "numMainQBins": "",
        "numWingQBins": "",
        "AnnularAngleBin": "1",
        "Qmin": "0.003",
        "Qmax": "",
        "useErrorWeighting": false,
        "useMaskBackTubes": false,
        "wavelength": "",
        "wavelengthSpread": "",
        "overlapStitchQmin": "0.075",
        "overlapStitchQmax": "0.095",
        "timeslice": false,
        "timesliceinterval": "200",
        "logslicename": "",
        "logslice": false,
        "logsliceinterval": "",
        "SampleToSi": "",
        "SampleDetectorDistance": ""
        }
    }"""

    # replace
    if sample_to_si_window_distance is not None:
        json_str = json_str.replace('"SampleToSi": ""',
                                    '"SampleToSi": "{}"'.format(sample_to_si_window_distance))

    if sample_to_detector_distance is not None:
        json_str = json_str.replace('"SampleDetectorDistance": ""',
                                    '"SampleDetectorDistance": "{}"'.format(sample_to_detector_distance))

    return json_str


def verify_reduction_results(sample_names, output_dir, gold_path, prefix):
    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_{}_reduction_log.hdf'.format(sample_name, prefix))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        compare_reduced_iq(output_log_file, gold_log_file)


def compare_reduced_iq(test_log_file, gold_log_file):
    """

    Parameters
    ----------
    test_log_file
    gold_log_file

    Returns
    -------

    """
    log_errors = list()

    for is_main_detector in [True, False]:
        vec_q_a, vec_i_a = get_iq1d(test_log_file, is_main=is_main_detector)
        vec_q_b, vec_i_b = get_iq1d(gold_log_file, is_main=is_main_detector)

        try:
            np.testing.assert_allclose(vec_q_a, vec_q_b)
            np.testing.assert_allclose(vec_i_a, vec_i_b)
            log_errors.append(None)
        except AssertionError as assert_err:
            log_errors.append(assert_err)
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
