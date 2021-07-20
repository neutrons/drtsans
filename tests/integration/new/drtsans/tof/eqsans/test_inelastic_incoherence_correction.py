import pytest
import os
from jsonschema.exceptions import ValidationError
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration, plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch
import numpy as np
import tempfile
from drtsans.dataobjects import save_i_of_q_to_h5, load_iq1d_from_h5
from typing import List, Any, Union, Dict
from matplotlib import pyplot as plt


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason='Required test data not available')
def test_parse_json():
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
        # TODO - expect to fail as elastic reference run 260159121 does not exist
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


# FIXME - better test method name
@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-25813/nexus/EQSANS_113915.nxs.h5'),
                    reason="Required test data not available")
def test_689_test2(reference_dir):
    from drtsans.tof.eqsans.correction_api import CorrectionConfiguration
    from drtsans.tof.eqsans.correction_api import ElasticReferenceRunSetup

    # Set up the configuration dict
    configuration = generate_configuration_with_correction()

    # Defaults and Validate
    input_config = reduction_parameters(configuration)

    # Check that inelastic incoherence config items were parsed
    # FIXME 777 and 785 - make this a method (centralized)
    if input_config['configuration'].get('fitInelasticIncoh'):
        select_min_incoh = input_config['configuration'].get('selectMinIncoh')
        correction = CorrectionConfiguration(True, select_min_incoh)
        elastic_run = input_config['configuration']['elasticReference'].get('runNumber')
        if elastic_run is not None:
            elastic_setup = ElasticReferenceRunSetup(elastic_run)
            correction.set_elastic_reference_run(elastic_setup)
    else:
        correction = None

    # Create temp output directory
    test_dir = tempfile.mkdtemp()
    print(f'[DEBUG 777] temp dir {test_dir} is created now')

    # 689 - this is only for debugging and TDD
    # base_name = ['baseline', 'newworkflow', 'correction'][2]
    # base_name = f'EQSANS_113915_{base_name}'
    base_name = 'EQSANS_113915_Incoh_1d'

    assert os.path.exists(test_dir), f'Output dir {test_dir} does not exit'
    configuration['configuration']['outputDir'] = test_dir
    configuration['outputFileName'] = base_name
    configuration['dataDirectories'] = test_dir

    # validate and clean configuration
    input_config = reduction_parameters(configuration)
    loaded = load_all_files(input_config)
    if base_name.endswith('baseline'):
        # Baseline
        correction = None
    else:
        pass

    # Reduce
    reduction_output = reduce_single_configuration(loaded, input_config,
                                                   incoherence_correction_setup=correction)

    if base_name.endswith('baseline'):
        # output reduced result as gold files
        export_reduction_output(reduction_output, output_dir=test_dir, prefix='EQSANS_113915')
    else:
        # Verify reduced workspace (nothing to do correction)
        reduced_data_nexus = os.path.join(test_dir, f'{base_name}.nxs')
        assert os.path.exists(reduced_data_nexus), f'Expected {reduced_data_nexus} does not exist'

        # verify with gold data
        # gold_dir = os.path.join(reference_dir.new.eqsans, 'baseline')
        gold_dir = reference_dir.new.eqsans
        gold_file = os.path.join(gold_dir, f'EQSANS_113915_baseline.nxs')
        verify_reduced_workspace(test_processed_nexus=reduced_data_nexus, gold_processed_nexus=gold_file,
                                 ws_prefix='new')


def verify_reduced_workspace(test_processed_nexus, gold_processed_nexus, ws_prefix):
    """Verify reduced result by verified expected result

    Parameters
    ----------
    test_processed_nexus: str
        NeXus file from test to verify
    gold_processed_nexus: str
        NeXus file containing the expected reduced result to verify against
    ws_prefix: str
        prefix for Mantid workspace that the

    """
    assert os.path.exists(gold_processed_nexus), f'Gold file {gold_processed_nexus} cannot be found'
    gold_ws = LoadNexusProcessed(Filename=gold_processed_nexus, OutputWorkspace=f'{ws_prefix}_gold')
    test_ws = LoadNexusProcessed(Filename=test_processed_nexus, OutputWorkspace=f'{ws_prefix}_test')
    r = CheckWorkspacesMatch(Workspace1=gold_ws, Workspace2=test_ws)
    if r != 'Success':
        assert gold_ws.getNumberHistograms() == test_ws.getNumberHistograms(), \
            f'Histograms: {gold_ws.getNumberHistograms()} != {test_ws.getNumberHistograms()}'
        assert gold_ws.readY(0).shape == test_ws.readY(0).shape, \
            f'Number of wavelength: {gold_ws.readY(0).shape} != {test_ws.readY(0).shape}'
        assert gold_ws.readX(0).shape == test_ws.readX(0).shape, \
            f'Histogram or point data: {gold_ws.readX(0).shape} != {test_ws.readX(0).shape}'
        gold_x_array = gold_ws.extractX()
        test_x_array = test_ws.extractX()
        assert gold_x_array.shape == test_x_array.shape
        np.testing.assert_allclose(gold_ws.extractX(), test_ws.extractX())
        np.testing.assert_allclose(gold_ws.extractY(), test_ws.extractY())
        np.testing.assert_allclose(gold_ws.extractE(), test_ws.extractE())


# FIXME - combine all the export_reduction_output() and add to other module
def export_reduction_output(reduction_output: List[Any], output_dir: Union[None, str] = None, prefix: str = ''):
    """Export the reduced I(Q) and I(Qx, Qy) to  hdf5 files
    """
    if output_dir is None:
        output_dir = os.getcwd()

    for section_index, section_output in enumerate(reduction_output):
        # 1D (list of IQmod)
        iq1ds = section_output.I1D_main
        for j_index, iq1d in enumerate(iq1ds):
            save_i_of_q_to_h5(iq1d, os.path.join(output_dir, f'{prefix}iq1d_{section_index}_{j_index}.h5'))
        # 2D (IQazimuthal)
        iq2d = section_output.I2D_main
        save_i_of_q_to_h5(iq2d, os.path.join(output_dir, f'{prefix}iq2d_{section_index}.h5'))


if __name__ == "__main__":
    pytest.main([__file__])
