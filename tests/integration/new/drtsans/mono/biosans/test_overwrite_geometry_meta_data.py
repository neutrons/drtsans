# Integration test for overwriting instrument geometry related meta data for BIO-SANS
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
    sample_thick = ['0.1', '0.1']  # both thickness are 0.1 cm
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
        reduction_input["thickness"] = sample_thick[i]  # thickness is in unit cm
        loaded = load_all_files(reduction_input)
        out = reduce_single_configuration(loaded, reduction_input)
        assert out

    end_time = time.time()
    print(end_time - start_time)


def skip_test_no_overwrite():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    # Set test and run
    json_file = generate_testing_json()
    output_dir = '/tmp/meta_overwrite_test1'

    reduce_biosans_data(json_file, output_dir)


def skip_test_over_write_both():
    """

    Returns
    -------

    """
    # Set test and output directory
    output_dir = '/tmp/bio_meta_overwrite_test1'

    # Json string
    input_json = generate_testing_json()
    # replace
    input_json = input_json.replace('"SampleToSi": ""', '"SampleToSi": "???.???"')
    input_json = input_json.replace('"SampleDetectorDistance": ""', '"SampleDetectorDistance": "??.??"')

    reduce_biosans_data(input_json, output_dir)


def generate_testing_json():
    """Generating testing JSON

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
        "sampleApertureSize":"14",
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
        "EvenDecades": false,
        "WedgeMinAngles": "-30, 60",
        "WedgeMaxAngles": "30, 120",
        "numMainQBins": "50",
        "numWingQBins": "50",
        "AnnularAngleBin": "1",
        "Qmin": "0.007",
        "Qmax": "0.85",
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

    return json_str


if __name__ == '__main__':
    pytest.main(__file__)
