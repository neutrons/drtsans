import pytest
import os
import json
from mantid.simpleapi import mtd
from drtsans.mono.biosans import load_all_files
from drtsans.mono.biosans import reduction_parameters
from drtsans.geometry import sample_detector_distance
from drtsans.samplelogs import SampleLogs


def test_load_all_files(reference_dir):
    """Standard reduction workflow

    Parameters
    ----------'
    reference_dir: str
        fixture path to testing data

    Returns
    -------

    """
    # Create JSON
    json_str = generate_test_json(os.path.join(reference_dir.new.biosans, 'overwrite_gold'))
    # Load JSON for configuration
    reduction_input = json.loads(json_str)

    # Set output and inputs
    output_dir = '/tmp/test_bio_load/'
    nexus_dir = reference_dir.new.biosans

    # Set up (testing) runs
    sample_names = ['csmb_ecoli1h_n2']
    samples = ['5709']
    samples_trans = samples
    backgrounds = ['5715']
    backgrounds_trans = backgrounds

    # chekcing if output directory exists, if it doesn't, creates the folder
    reduction_input["configuration"]["outputDir"] = output_dir

    reduction_input["sample"]["runNumber"] = samples[0]
    reduction_input["sample"]["transmission"]["runNumber"] = samples_trans[0]
    reduction_input["background"]["runNumber"] = backgrounds[0]
    reduction_input["background"]["transmission"]["runNumber"] = backgrounds_trans[0]
    reduction_input["outputFileName"] = sample_names[0]
    reduction_input["dataDirectories"] = nexus_dir

    # convert
    reduction_input = reduction_parameters(reduction_input, validate=True)
    load_all_files(reduction_input, path=nexus_dir, prefix='BioTestLoadAll')

    beam_center_run = mtd['BioTestLoadAll_BIOSANS_1322_raw_histo']
    dark_run = mtd['BioTestLoadAll_BIOSANS_1383_raw_histo']
    empty_trans_run = mtd['BioTestLoadAll_BIOSANS_5705_raw_histo']
    sample_run = mtd['BioTestLoadAll_BIOSANS_5709_raw_histo']
    bkgd_run = mtd['BioTestLoadAll_BIOSANS_5715_raw_histo']

    # Verify sample to si-window distance by checking the sample position with default setup
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, beam_center_run, bkgd_run, empty_trans_run]:
        sample_pos_z = ws.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(-0.12952, 0.000004), '{} has a wrong sample Si-window distance {}' \
                                                                  ''.format(str(ws), sample_pos_z)

    for ws in [dark_run]:
        sample_pos_z = ws.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(0.0000, 0.000004), '{} has a wrong sample Si-window distance {}' \
                                                                ''.format(str(ws), sample_pos_z)

    # Verify sample to detector distance with default setup:
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/542#note_156296
    for ws in [sample_run, beam_center_run, bkgd_run, empty_trans_run]:
        # reset SDD with sample run
        sdd_value = sample_detector_distance(ws, unit='m', search_logs=False)
        assert sdd_value == pytest.approx(14.31, 0.004), '{} has a wrong SDD {}'.format(str(ws), sdd_value)

    for ws in [dark_run]:
        # EPICS recorded SDD with sample run
        sdd_value = sample_detector_distance(ws, unit='m', search_logs=False)
        assert sdd_value == pytest.approx(7.00, 1E-6), '{} has a wrong SDD {}'.format(str(ws), sdd_value)

    # Verify smearing pixel size x and smearing pixel size y
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/541: default
    for ws_index, ws in enumerate([sample_run, beam_center_run, bkgd_run, empty_trans_run, dark_run]):
        sample_log_i = SampleLogs(ws)
        pixel_size_x = sample_log_i['smearingPixelSizeX'].value
        pixel_size_y = sample_log_i['smearingPixelSizeY'].value
        assert pixel_size_x == pytest.approx(1.23456 * 1.E-3, 1.E-7), \
            '{}-th workspace: Pixel size X {} (m) shall be equal to 1.23456 mm'.format(ws_index, pixel_size_x)
        assert pixel_size_y == pytest.approx(2.34567 * 1.E-3, 1.E-7), \
            '{}-th workspace: Pixel size X {} (m) shall be equal to 2.34567 mm'.format(ws_index, pixel_size_y)


def generate_test_json(sens_nxs_dir):
    """

    Parameters
    ----------
    sens_nxs_dir: str
        path to sensitivity nexus file

    Returns
    -------
    dict
        configuration JSON

    """
    specs = """{
     "instrumentName": "BIOSANS",
     "iptsNumber": "24740",
     "sample": {
        "runNumber": "4822",
        "thickness": "0.1",
        "transmission": {
            "runNumber": "4822"
        }
     },
     "outputFileName": "CG3_4822",
     "background": {
         "runNumber": "4821",
         "transmission": {
             "runNumber": "4821"
         }
     },
     "beamCenter": {
         "runNumber": "1322"
     },
     "emptyTransmission": {
         "runNumber": "5705"
     },
     "configuration": {
         "outputDir": "/HFIR/CG3/shared/UserAcceptance/override/test1",
         "sampleApertureSize": "14",
         "sourceApertureDiameter": "",
         "maskFileName": "",
         "useDefaultMask": true,
         "defaultMask": ["{'Pixel':'1-12,244-256'}", "{'Bank':'21-24,45-48'}"],
         "blockedBeamRunNumber": "",
         "darkMainFileName": "CG3_1383.nxs",
         "darkWingFileName": "CG3_1383.nxs",
         "sensitivityMainFileName": "AAA",
         "sensitivityWingFileName": "BBB",
         "absoluteScaleMethod": "standard",
         "DBScalingBeamRadius": "",
         "StandardAbsoluteScale": "0.0055e-8",
         "normalization": "Monitor",
         "sampleOffset": "",
         "useSolidAngleCorrection": true,
         "useThetaDepTransCorrection": true,
         "mmRadiusForTransmission": "",
         "numMainQxQyBins": "100",
         "numWingQxQyBins": "100",
         "1DQbinType": "scalar",
         "QbinType": "log",
         "useLogQBinsEvenDecade": false,
         "LogQBinsPerDecadeMain": 20,
         "LogQBinsPerDecadeWing": 25,
         "WedgeMinAngles": "-30, 60",
         "WedgeMaxAngles": "30, 120",
         "numMainQBins": "",
         "numWingQBins": "",
         "AnnularAngleBin": "1",
         "QminMain": 0.003,
         "QminWing": 0.003,
         "QmaxMain": "",
         "QmaxWing": "",
         "useErrorWeighting": false,
         "useMaskBackTubes": false,
         "wavelength": "",
         "wavelengthSpread": "",
         "overlapStitchQmin": "0.075",
         "overlapStitchQmax": "0.095",
         "useTimeSlice": false,
         "timeSliceInterval": 200,
         "logSliceName": "",
         "useLogSlice": false,
         "logSliceInterval": "",
         "sampleToSi": "200.52",
         "sampleDetectorDistance": "14.31",
         "smearingPixelSizeX": "1.23456",
         "smearingPixelSizeY": "2.34567"
         }
     }"""

    # replace values
    main_sens = os.path.join(sens_nxs_dir, 'sens_f4829m7p0_TDC_SAC.h5')
    specs = specs.replace('AAA', main_sens)
    wing_sens = os.path.join(sens_nxs_dir, 'sens_f4835w3p2_TDC_SAC.h5')
    specs = specs.replace('BBB', wing_sens)

    return specs


if __name__ == '__main__':
    pytest.main([__file__])