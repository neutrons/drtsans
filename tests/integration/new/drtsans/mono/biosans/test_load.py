import pytest
import os
import json
from mantid.simpleapi import mtd
from drtsans.mono.biosans import load_all_files
# from drtsans.geometry import sample_detector_distance


def test_load_all_files(reference_dir):
    """Standard reduction workflow

    Parameters
    ----------'
    reference_dir: str
        fixture path to testing data

    Returns
    -------

    """
    # Set output and inputs
    output_dir = '/tmp/test_bio_load/'
    nexus_dir = reference_dir.new.biosans

    # Create JSON
    json_str = generate_test_json(os.path.join(reference_dir.new.biosans, 'overwrite_gold'))

    # Set up (testing) runs
    sample_names = ['csmb_ecoli1h_n2']
    samples = ['5709']
    samples_trans = samples
    backgrounds = ['5715']
    backgrounds_trans = backgrounds

    # Load JSON for configuration
    reduction_input = json.loads(json_str)

    # chekcing if output directory exists, if it doesn't, creates the folder
    reduction_input["configuration"]["outputDir"] = output_dir

    reduction_input["runNumber"] = samples[0]
    reduction_input["transmission"]["runNumber"] = samples_trans[0]
    reduction_input["background"]["runNumber"] = backgrounds[0]
    reduction_input["background"]["transmission"]["runNumber"] = backgrounds_trans[0]
    reduction_input["outputFilename"] = sample_names[0]
    load_all_files(reduction_input, path=nexus_dir, prefix='BioTestLoadAll')

    # Outputs
    output_ws_names = mtd.getObjectNames()
    for ws_name in output_ws_names:
        print('DDD  {}'.format(ws_name))

    assert 1 == 5


def generate_test_json(sens_nxs_dir):
    """

    Parameters
    ----------
    sens_nxs_dir: str
        path to sensitivity nexus file

    Returns
    -------
    str
        configuration JSON

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

    # replace values
    #  '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4829m7p0_TDC_SAC.h5'
    main_sens = os.path.join(sens_nxs_dir, 'sens_f4829m7p0_TDC_SAC.h5')
    json_str = json_str.replace('/HFIR/CG3/shared/Cycle486/sens_f4829m7p0_TDC_SAC.h5', main_sens)

    # '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4835w3p2_TDC_SAC.h5'
    wing_sens = os.path.join(sens_nxs_dir, 'sens_f4835w3p2_TDC_SAC.h5')
    json_str = json_str.replace('/HFIR/CG3/shared/Cycle486/sens_f4835w3p2_TDC_SAC.h5', wing_sens)

    return json_str


if __name__ == '__main__':
    pytest.main([__file__])
