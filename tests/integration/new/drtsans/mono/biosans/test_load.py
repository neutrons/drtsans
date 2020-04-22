import pytest
import os
import json
from mantid.simpleapi import mtd
from drtsans.mono.biosans import load_all_files
from drtsans.geometry import sample_detector_distance


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

    beam_center_run = mtd['BioTestLoadAll_CG3_1322_raw_histo']
    dark_main_run = mtd['BioTestLoadAll_CG3_1383_raw_histo']
    dark_wing_run = mtd['BioTestLoadAll_CG3_5705_raw_histo']
    sample_run = mtd['BioTestLoadAll_CG3_5709_raw_histo']
    bkgd_run = mtd['BioTestLoadAll_CG3_5715_raw_histo']

    # Verify sample to si-window distance by checking the sample position with default setup
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, beam_center_run, bkgd_run]:
        sample_pos_z = ws.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(-0.12952, 0.000004), '{} has a wrong sample Si-window distance {}' \
                                                                  ''.format(str(ws), sample_pos_z)

    for ws in [dark_main_run, dark_wing_run]:
        sample_pos_z = ws.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(-0.071, 0.000004), '{} has a wrong sample Si-window distance {}' \
                                                                ''.format(str(ws), sample_pos_z)

    # Verify sample to detector distance with default setup:
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, beam_center_run, bkgd_run]:
        # reset SDD with sample run
        sdd_value = sample_detector_distance(ws, unit='m', search_logs=False)
        assert sdd_value == pytest.approx(14.31, 0.004), '{} has a wrong SDD {}'.format(str(ws), sdd_value)

    for ws in [dark_wing_run, dark_main_run]:
        # EPICS recorded SDD with sample run
        sdd_value = sample_detector_distance(ws, unit='m', search_logs=False)
        assert sdd_value == pytest.approx(10.00, 0.004), '{} has a wrong SDD {}'.format(str(ws), sdd_value)


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
         "SampleToSi": "200.52",
         "SampleDetectorDistance": "14.31"
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
