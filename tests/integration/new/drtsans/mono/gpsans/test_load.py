import pytest
import os
import json
from mantid.simpleapi import mtd
from drtsans.mono.gpsans import load_all_files
from drtsans.geometry import sample_detector_distance
# from drtsans.mono.gpsans import load_and_split


def test_load_all_files(reference_dir):
    """Standard reduction workflow

    Parameters
    ----------'
    data_dir
    json_file
    output_dir

    Returns
    -------

    """
    # Set output
    output_dir = '/tmp/test_load/'

    # Create JSON file
    json_file = generate_test_json()

    # USER Input here with scan numbers etc.
    samples = ['9166']
    samples_trans = ['9178']
    sample_thick = ['0.1']
    bkgd = ['9165']
    bkgd_trans = ['9177']

    # Sample names for output
    sample_names = ["Al4"]

    # Import JSON
    reduction_input = json.loads(json_file)

    # set output directory
    reduction_input["configuration"]["outputDir"] = output_dir

    reduction_input["runNumber"] = samples[0]
    reduction_input["transmission"]["runNumber"] = samples_trans[0]
    reduction_input["background"]["runNumber"] = bkgd[0]
    reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[0]
    reduction_input["outputFilename"] = sample_names[0]
    reduction_input["thickness"] = sample_thick[0]
    load_all_files(reduction_input, prefix='GP_TEST_LOAD', load_params=None, path=reference_dir.new.gpsans)

    sample_run = mtd['GP_TEST_LOAD_CG2_9166_raw_histo']
    sample_trans_run = mtd['GP_TEST_LOAD_CG2_9178_raw_histo']
    bkgd_run = mtd['GP_TEST_LOAD_CG2_9165_raw_histo']
    bkgd_trans_run = mtd['GP_TEST_LOAD_CG2_9177_raw_histo']

    # Verify sample to si-window distance by checking the sample position with default setup
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, sample_trans_run, bkgd_run, bkgd_trans_run]:
        sample_pos_z = ws.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(0.23456, 0.000004), '{} has a wrong sample Si-window distance {}' \
                                                                 ''.format(str(ws), sample_pos_z)

    # Verify sample to detector distance with default setup:
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, sample_trans_run, bkgd_run, bkgd_trans_run]:
        # reset SDD with sample run
        sdd_value = sample_detector_distance(ws, unit='m', search_logs=False)
        assert sdd_value == pytest.approx(32.11, 0.004), '{} has a wrong SDD {}'.format(str(ws), sdd_value)


def generate_test_json():
    """Generate testing JSON

    Returns
    -------
    str
        JSON string

    """
    json_str = """
    {    "instrumentName": "CG2",
    "iptsNumber": "24727",
    "runNumber": "",
    "thickness": "0.1",
    "outputFilename": "",
    "transmission": {
        "runNumber": "",
        "value": ""
    },
    "background": {
        "runNumber": "",
        "transmission": {
            "runNumber": "",
            "value": ""
        }
    },
    "beamCenter": {
        "runNumber": "9177"
    },
    "emptyTrans": {
        "runNumber": "9177"
    },
    "configuration": {
        "outputDir": "/HFIR/CG2/shared/UserAcceptance/overwrite_meta/test1/",
        "sampleApertureSize": "",
        "sourceApertureDiameter": "",
        "maskFileName": "",
        "useMaskFileName": false,
        "useDefaultMask": true,
        "DefaultMask":["{'Pixel':'1-10,247-256'}"],
        "useBlockedBeam": false,
        "BlockBeamFileName":"",
        "useDarkFileName": false,
        "darkFileName": "",
        "useSensitivityFileName": true,
        "sensitivityFileName": "SensFileToReplace",
        "UseBarScan": false,
        "BarScanFileName": "",
        "absoluteScaleMethod": "direct_beam",
        "DBScalingBeamRadius": "40",
        "StandardAbsoluteScale": "1",
        "normalization": "Monitor",
        "sampleOffset": "",
        "useSampleOffset": false,
        "useSolidAngleCorrection": true,
        "useThetaDepTransCorrection": true,
        "mmRadiusForTransmission": "40",
        "numQxQyBins": "150",
        "1DQbinType": "scalar",
        "QbinType": "log",
        "LogQBinsPerDecade": "33",
        "LogQBinsDecadeCenter": false,
        "LogQBinsEvenDecade": true,
        "WedgeMinAngles": "-30, 60",
        "WedgeMaxAngles": "30, 120",
        "numQBins": "",
        "AnnularAngleBin": "",
        "Qmin": "",
        "Qmax": "",
        "useErrorWeighting": false,
        "useMaskBackTubes": false,
        "wavelength": "",
        "wavelengthSpread": "",
        "timeslice": false,
        "timesliceinterval": "",
        "logslicename": "",
        "logslice": false,
        "logsliceinterval": "",
        "SampleToSi": "234.56",
        "SampleDetectorDistance": "32.11"
    }}
    """

    # replace for the correct sensitivity file
    sens_nxs_dir = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/gpsans/'
    sens_path = os.path.join(sens_nxs_dir, 'sens_c486_noBar.nxs')
    json_str = json_str.replace('SensFileToReplace', sens_path)

    return json_str


if __name__ == '__main__':
    pytest.main([__file__])
