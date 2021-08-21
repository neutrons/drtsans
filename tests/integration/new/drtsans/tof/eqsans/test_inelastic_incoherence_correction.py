import pytest
import os
from jsonschema.exceptions import ValidationError
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration)  # noqa E402
from drtsans.dataobjects import _Testing
import json
import tempfile
from typing import Tuple, Dict
from drtsans.dataobjects import load_iq1d_from_h5, load_iq2d_from_h5

@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason='Required test data not available')
def test_parse_json():
    """Test the JSON to dictionary
    """
    elastic_reference_run = "124680"
    elastic_reference_bkgd_run = "null"
    # Specify JSON input
    reduction_input = {
        "instrumentName": "EQSANS",
        "iptsNumber": "26015",
        "sample": {
            "runNumber": "115363",
            "thickness": "1.0"
        },
        "background": {
            "runNumber": "115361",
            "transmission": {
                "runNumber": "115357",
                "value": ""
            }
        },
        "emptyTransmission": {
            "runNumber": "115356",
            "value": ""
        },
        "beamCenter": {
            "runNumber": "115356"
        },
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "outputDir": "/path/to/nowhere",
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": ("/SNS/EQSANS/"
                                     "IPTS-24769/shared/EQSANS_110943.out"),
            "sensitivityFileName": ("/SNS/EQSANS/"
                                    "shared/NeXusFiles/EQSANS/"
                                    "2020A_mp/Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs"),
            "numQBins": "100",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": True,
            "fitInelasticIncoh": True,
            "elasticReference": {
                "runNumber": elastic_reference_run,
                "thickness": "1.0",
                "transmission": {
                    "runNumber": "null",
                    "value": "1.0"
                }
            },
            "elasticReferenceBkgd": {
                "runNumber": elastic_reference_bkgd_run,
                "transmission": {
                    "runNumber": "null",
                    "value": "0.9"
                }
            },
            "selectMinIncoh": True
        }
    }

    # Validate
    input_config = reduction_parameters(reduction_input)

    # Check that inelastic incoherence config items were parsed
    assert input_config['configuration'].get('fitInelasticIncoh')
    assert input_config['configuration']['elasticReference'].get('runNumber') == elastic_reference_run
    assert input_config['configuration'].get('selectMinIncoh')


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason='Required test data not available')
def test_parse_invalid_json():
    """Test the JSON to dictionary
    """
    invalid_run_num = "260159121"
    valid_run_num = "115363"
    # Specify JSON input
    reduction_input = {
        "instrumentName": "EQSANS",
        "iptsNumber": "26015",
        "sample": {
            "runNumber": "115363",
            "thickness": "1.0"
        },
        "background": {
            "runNumber": "115361",
            "transmission": {
                "runNumber": "115357",
                "value": ""
            }
        },
        "emptyTransmission": {
            "runNumber": "115356",
            "value": ""
        },
        "beamCenter": {
            "runNumber": "115356"
        },
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "outputDir": "/path/to/nowhere",
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": ("/SNS/EQSANS/"
                                     "IPTS-24769/shared/EQSANS_110943.out"),
            "sensitivityFileName": ("/SNS/EQSANS/"
                                    "shared/NeXusFiles/EQSANS/"
                                    "2020A_mp/Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs"),
            "numQBins": "100",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": True,
            "fitInelasticIncoh": True,
            "elasticReference": {
                "runNumber": invalid_run_num,
                "thickness": "1.0",
                "transmission": {
                    "runNumber": valid_run_num,
                    "value": "0.9"
                }
            },
            "elasticReferenceBkgd": {
                "runNumber": valid_run_num,
                "transmission": {
                    "runNumber": valid_run_num,
                    "value": "0.9"
                }
            },
            "selectMinIncoh": True
        }
    }

    # Validate
    with pytest.raises(ValidationError):
        # expect to fail as elastic reference run 260159121 does not exist
        reduction_parameters(reduction_input)

    # Respecify to use a valid run
    # json_str.replace('260159121', '26015')
    reduction_input['configuration']['elasticReference']['runNumber'] = valid_run_num
    # Defaults and Validate
    input_config = reduction_parameters(reduction_input)

    # Check that inelastic incoherence config items were parsed
    assert input_config['configuration'].get('fitInelasticIncoh')
    assert input_config['configuration']['elasticReference'].get('runNumber') == valid_run_num
    assert input_config['configuration'].get('selectMinIncoh')


def generate_configuration_with_correction(output_dir: str = '/tmp/') -> Dict:
    """Generate configuration dictionary (JSON) from test 2 in issue 689

    Source: https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689  *Test 2*

    Returns
    -------
    ~dict
        Reduction donfiguration JSON

    """
    reduction_configuration = {
        "schemaStamp": "2020-05-21T17:18:11.528041",
        "instrumentName": "EQSANS",
        "iptsNumber": "24876",
        "sample": {
            "runNumber": "113915",
            "thickness": 1,
            "transmission": {"runNumber": "113914", "value": ""}},
        "background": {
            "runNumber": "113919",
            "transmission": {"runNumber": "113918", "value": ""}},
        "emptyTransmission": {"runNumber": "113682", "value": None},
        "beamCenter": {"runNumber": "113682"},
        "outputFileName": "water65D_2o5m2o5a_full",
        "configuration": {
            "outputDir": f'{output_dir}',
            "useTimeSlice": False,
            "timeSliceInterval": 300,
            "useLogSlice": False,
            "logSliceName": None,
            "logSliceInterval": 10,
            "cutTOFmin": 500,
            "cutTOFmax": 2000,
            "wavelengthStep": 0.1,
            "wavelengthStepType": "constant Delta lambda",
            "sampleOffset": 314.5,
            "useDetectorOffset": True,
            "detectorOffset": 80,
            "sampleApertureSize": 10,
            "sourceApertureDiameter": None,
            "usePixelCalibration": None,
            "maskFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2020A_mp/mask_4m_extended.nxs",
            "useDefaultMask": True,
            "defaultMask": None,
            "useMaskBackTubes": False,
            "darkFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2020A_mp/EQSANS_113569.nxs.h5",
            "normalization": "Total charge",
            "fluxMonitorRatioFile": None,
            "beamFluxFileName": "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample",
            "sensitivityFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2020A_mp/"
                                   "Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs",
            "useSolidAngleCorrection": True,
            "useThetaDepTransCorrection": True,
            "mmRadiusForTransmission": 25,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 1.0,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 100,
            "LogQBinsPerDecade": None,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "autoWedgeQmin": None,
            "autoWedgeQmax": None,
            "autoWedgeQdelta": None,
            "autoWedgeAzimuthalDelta": None,
            "autoWedgePeakWidth": None,
            "autoWedgeBackgroundWidth": None,
            "autoWedgeSignalToNoiseMin": None,
            "AnnularAngleBin": 5.0,
            "Qmin": None,
            "Qmax": None,
            "useErrorWeighting": True,
            "smearingPixelSizeX": None,
            "smearingPixelSizeY": None,
            "useSubpixels": None,
            "subpixelsX": None,
            "subpixelsY": None,
            "useSliceIDxAsSuffix": None,
            # inelastic/incoherent correction
            "fitInelasticIncoh": True,
            "elasticReference": {
                "runNumber": None,
                "thickness": "1.0",
                "transmission": {
                    "runNumber": None,
                    "value": "0.9"
                }
            },
            "elasticReferenceBkgd": {
                "runNumber": None,
                "transmission": {
                    "runNumber": None,
                    "value": "0.9"
                }
            },
            "selectMinIncoh": True
        }
    }

    return reduction_configuration


@pytest.mark.skipif(not os.path.exists('/SNS/users/pf9/etc/'),
                    reason="Test is too long for build server")
def test_incoherence_correction(reference_dir):
    """Test incoherence correction without elastic correction
    """
    # Set up the configuration dict
    configuration = generate_configuration_with_correction()

    # Create temp output directory
    test_dir = tempfile.mkdtemp()
    base_name = 'EQSANS_113915_Incoh_1d'

    assert os.path.exists(test_dir), f'Output dir {test_dir} does not exit'
    configuration['configuration']['outputDir'] = test_dir
    configuration['outputFileName'] = base_name
    configuration['dataDirectories'] = test_dir

    # validate and clean configuration
    input_config = reduction_parameters(configuration)
    loaded = load_all_files(input_config)

    # Reduce
    reduction_output = reduce_single_configuration(loaded, input_config,
                                                   not_apply_incoherence_correction=False)

    # Gold data directory
    gold_dir = os.path.join(reference_dir.new.eqsans, 'gold_data/Incoherence_Corrected_113915/')
    assert os.path.exists(gold_dir), f'Gold/expected data directory {gold_dir} does not exist'

    # Verify with gold data
    gold_file_dict = dict()
    for frame_index in range(1):
        iq1d_h5_name = os.path.join(gold_dir, f'EQSANS_11395iq1d_{frame_index}_0.h5')
        gold_file_dict[1, frame_index, 0] = iq1d_h5_name
        iq2d_h5_name = os.path.join(gold_dir, f'EQSANS_11395iq2d_{frame_index}.h5')
        gold_file_dict[2, frame_index] = iq2d_h5_name
        assert os.path.exists(iq1d_h5_name) and os.path.exists(iq2d_h5_name), f'{iq1d_h5_name} and/or {iq2d_h5_name}' \
                                                                              f'do not exist'

    # Verify
    verify_binned_iq(gold_file_dict, reduction_output)


@pytest.mark.skipif(not os.path.exists('/SNS/users/pf9/etc/'),
                    reason="Test is too long for build server")
def test_incoherence_correction_elastic_normalization(reference_dir):
    """Test incoherence correction with elastic correction
    """
    # Set up the configuration dict
    config_json_file = os.path.join(reference_dir.new.eqsans, 'test_incoherence_correction/agbe_125707_test1.json')
    assert os.path.exists(config_json_file), f'Test JSON file {config_json_file} does not exist.'
    with open(config_json_file, 'r') as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Create temp output directory
    test_dir = tempfile.mkdtemp()
    base_name = 'EQSANS_113915_Incoh_1d'

    assert os.path.exists(test_dir), f'Output dir {test_dir} does not exit'
    configuration['configuration']['outputDir'] = test_dir
    configuration['outputFileName'] = base_name
    configuration['dataDirectories'] = test_dir

    # validate and clean configuration
    input_config = reduction_parameters(configuration)
    loaded = load_all_files(input_config)

    # Reduce
    reduction_output = reduce_single_configuration(loaded, input_config,
                                                   not_apply_incoherence_correction=False)

    # Gold data directory
    gold_dir = os.path.join(reference_dir.new.eqsans, 'gold_data/Incoherence_Corrected_113915/')
    assert os.path.exists(gold_dir), f'Gold/expected data directory {gold_dir} does not exist'

    # Verify with gold data
    gold_file_dict = dict()
    for frame_index in range(1):
        iq1d_h5_name = os.path.join(gold_dir, f'EQSANS_11395iq1d_{frame_index}_0.h5')
        gold_file_dict[1, frame_index, 0] = iq1d_h5_name
        iq2d_h5_name = os.path.join(gold_dir, f'EQSANS_11395iq2d_{frame_index}.h5')
        gold_file_dict[2, frame_index] = iq2d_h5_name
        assert os.path.exists(iq1d_h5_name) and os.path.exists(iq2d_h5_name), f'{iq1d_h5_name} and/or {iq2d_h5_name}' \
                                                                              f'do not exist'

    # Verify
    verify_binned_iq(gold_file_dict, reduction_output)


def verify_binned_iq(gold_file_dict: Dict[Tuple, str], reduction_output):
    """Verify reduced I(Q1D) and I(qx, qy) by expected/gold data

    Parameters
    ----------
    gold_file_dict: ~dict
        dictionary for gold files
    reduction_output: ~list
        list of binned I(Q1D) and I(qx, qy)

    """
    num_frames_gold = len(gold_file_dict) // 2
    assert num_frames_gold == len(reduction_output),\
        f'Frame numbers are different: gold = {len(gold_file_dict) // 2}; test = {len(reduction_output)}'

    for frame_index in range(num_frames_gold):
        # 1D
        iq1d_h5_name = gold_file_dict[1, frame_index, 0]
        gold_iq1d = load_iq1d_from_h5(iq1d_h5_name)
        _Testing.assert_allclose(reduction_output[frame_index].I1D_main[0], gold_iq1d)

        # 2D
        iq2d_h5_name = gold_file_dict[2, frame_index]
        gold_iq2d = load_iq2d_from_h5(iq2d_h5_name)
        _Testing.assert_allclose(reduction_output[frame_index].I2D_main, gold_iq2d)


if __name__ == "__main__":
    pytest.main([__file__])
