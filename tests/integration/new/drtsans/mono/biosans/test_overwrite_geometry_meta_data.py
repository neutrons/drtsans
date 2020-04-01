# Integration test for overwriting instrument geometry related meta data for BIO-SANS
# From round 3 4822
# Tests are
import pytest
import json
import os
from drtsans.mono.biosans import load_all_files, reduce_single_configuration
import time
from mantid.simpleapi import mtd


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
    sample_names = ['csmb_ecoli1h_org', 'insect1hTime_org']
    samples = ['5709', '5712']
    samples_trans = samples
    backgrounds = ['5715', '5715']
    backgrounds_trans = backgrounds

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
        loaded = load_all_files(reduction_input)
        out = reduce_single_configuration(loaded, reduction_input)
        assert out

    end_time = time.time()
    print(end_time - start_time)


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


def skip_test_over_write_both():
    """

    Returns
    -------

    """
    # Set test and output directory
    json_file = generate_testing_json(None, None)
    output_dir = '/tmp/bio_meta_overwrite_test4'

    # run
    reduce_biosans_data(json_file, output_dir)


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


if __name__ == '__main__':
    pytest.main(__file__)
