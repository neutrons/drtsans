import pytest
import os
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch
import numpy as np
import tempfile


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason='Required test data not available')
def failed_test_parse_json():
    """Test the JSON to dictionary
    """
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
            "elasticReference": "260159121",
            "selectMinIncoh": True
        }
    }

    # Validate
    with pytest.raises(RuntimeError):
        # TODO - expect to fail as elastic reference run 260159121 does not exist
        reduction_parameters(reduction_input)

    # Respecify to use a valid run
    # json_str.replace('260159121', '26015')
    reduction_input['configuration']['elasticReference'] = "115363"
    # Defaults and Validate
    input_config = reduction_parameters(reduction_input)

    # Check that inelastic incoherence config items were parsed
    assert input_config['configuration'].get('fitInelasticIncoh')
    assert input_config['configuration'].get('elasticReference') == '115363'
    assert input_config['configuration'].get('selectMinIncoh')


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason="Required test data not available")
def failed_test_correct_without_elastic(reference_dir):

    # Set up the configuration dict
    configuration = {
        'instrumentName': 'EQSANS',
        'iptsNumber': 26015,
        'sample': {
            'runNumber': 115363,
            'thickness': 1.0
        },
        'background': {
            'runNumber': 115361,
            'transmission': {
                'runNumber': 115357,
                'value': ''
            }
        },
        'emptyTransmission': {
            'runNumber': 115356,
            'value': ''
        },
        'beamCenter': {
            'runNumber': 115356
        },
        'outputFileName': 'test_wavelength_step',
        'configuration': {
            'cutTOFmax': '1500',
            'wavelengthStepType': 'constant Delta lambda/lambda',
            'sampleApertureSize': '10',
            'fluxMonitorRatioFile': ('/SNS/EQSANS/'
                                     'IPTS-24769/shared/EQSANS_110943.out'),
            'sensitivityFileName': ('/SNS/EQSANS/shared/NeXusFiles/'
                                    'EQSANS/2020A_mp/'
                                    'Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs'),
            'numQBins': 100,
            'WedgeMinAngles': '-30, 60',
            'WedgeMaxAngles': '30, 120',
            'AnnularAngleBin': '5',
            'useSliceIDxAsSuffix': True,
            # inelastic/incoherent correction
            "fitInelasticIncoh": True,
            "elasticReference": None,
            "selectMinIncoh": False,
        }
    }

    gold_file_dir = os.path.join(reference_dir.new.eqsans, '')
    assert os.path.exists(gold_file_dir), f'EQSANS test directory: {gold_file_dir} does not exist'
    gold_file_dir = os.path.join(reference_dir.new.eqsans, 'test_reduce/gold_data')
    assert os.path.exists(gold_file_dir), f'EQSANS gold data directory: {gold_file_dir} does not exist'

    # Create temp output directory
    # FIXME the following with block is skipped for quick test on the next with-block
    # with tempfile.TemporaryDirectory() as test_dir:
    #     configuration['configuration']['outputDir'] = test_dir
    #     configuration['outputFileName'] = 'test_wavelength_step_reg'
    #     configuration['dataDirectories'] = test_dir
    #     # validate and clean configuration
    #     input_config = reduction_parameters(configuration)
    #     loaded = load_all_files(input_config)
    #     reduced_value = reduce_single_configuration(loaded, input_config)

    # print(f'reduced: {reduced_value}')
    # print(f'reduced type = {type(reduced_value)}')
    # print(f'list size = {len(reduced_value)}')
    # for val in reduced_value:
    #     print(f'type = {type(val)}')

    from drtsans.tof.eqsans.correction_api import CorrectionConfiguration

    # Create output directory
    with tempfile.TemporaryDirectory() as test_dir:
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_reg'
        configuration['dataDirectories'] = test_dir

        # Test correction setup
        test_setup = CorrectionConfiguration(True, False)
        test_setup.debug_no_correction = True

        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config, incoherence_correction_setup=test_setup)


def verify_reduction(test_file, gold_file, ws_prefix):
    """Verify reduced result by verified expected result

    Parameters
    ----------
    test_file: str
        NeXus file from test to verify
    gold_file: str
        NeXus file containing the expected reduced result to verify against
    ws_prefix: str
        prefix for Mantid workspace that the

    """
    assert os.path.exists(gold_file), f'Gold file {gold_file} cannot be found'
    gold_ws = LoadNexusProcessed(Filename=gold_file, OutputWorkspace=f'{ws_prefix}_gold')
    test_ws = LoadNexusProcessed(Filename=test_file, OutputWorkspace=f'{ws_prefix}_test')
    r = CheckWorkspacesMatch(Workspace1=gold_ws, Workspace2=test_ws)
    if r != 'Success':
        assert gold_ws.getNumberHistograms() == test_ws.getNumberHistograms(),\
            f'Histograms: {gold_ws.getNumberHistograms()} != {test_ws.getNumberHistograms()}'
        assert gold_ws.readY(0).shape == test_ws.readY(0).shape,\
            f'Number of wavelength: {gold_ws.readY(0).shape} != {test_ws.readY(0).shape}'
        assert gold_ws.readX(0).shape == test_ws.readX(0).shape,\
            f'Histogram or point data: {gold_ws.readX(0).shape} != {test_ws.readX(0).shape}'
        gold_x_array = gold_ws.extractX()
        test_x_array = test_ws.extractX()
        assert gold_x_array.shape == test_x_array.shape
        np.testing.assert_allclose(gold_ws.extractX(), test_ws.extractX())
        np.testing.assert_allclose(gold_ws.extractY(), test_ws.extractY())
        np.testing.assert_allclose(gold_ws.extractE(), test_ws.extractE())


if __name__ == "__main__":
    pytest.main([__file__])
