import pytest
import os
import tempfile
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch
import numpy as np


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason="Required test data not available")
def test_wavelength_step(reference_dir):

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
            'useSliceIDxAsSuffix': True
        }
    }

    gold_file_dir = os.path.join(reference_dir.new.eqsans, '')
    assert os.path.exists(gold_file_dir), f'EQSANS test directory: {gold_file_dir} does not exist'
    gold_file_dir = os.path.join(reference_dir.new.eqsans, 'test_reduce/gold_data')
    assert os.path.exists(gold_file_dir), f'EQSANS gold data directory: {gold_file_dir} does not exist'

    # Create output directory
    with tempfile.TemporaryDirectory() as test_dir:
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_reg'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config)
        output_file_name = os.path.join(test_dir, 'test_wavelength_step_reg.nxs')
        assert os.path.isfile(output_file_name), f'Expected output file {output_file_name} does not exists'

        # verify_reduced_data
        gold_file = os.path.join(gold_file_dir, 'expected_wavelength_step_reg.nxs')
        exp_file = output_file_name
        verify_reduction(exp_file, gold_file, 'reg')

    with tempfile.TemporaryDirectory() as test_dir:
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_com'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        input_config['beamCenter']['method'] = 'center_of_mass'
        input_config['beamCenter']['com_centering_options'] = {'CenterX': 0., 'CenterY': 0., 'Tolerance': 0.00125}
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config)
        output_file_name = os.path.join(f'{test_dir}', 'test_wavelength_step_com.nxs')
        assert os.path.isfile(output_file_name),  f'Expected output file {output_file_name} does not exists'

        # verify_reduced_data
        gold_file = os.path.join(gold_file_dir, 'expected_wavelength_step_com.nxs')
        exp_file = output_file_name
        verify_reduction(exp_file, gold_file, 'com')

    with tempfile.TemporaryDirectory() as test_dir:
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_gauss'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        input_config['beamCenter']['method'] = 'gaussian'
        input_config['beamCenter']['gaussian_centering_options'] = {'theta': {'value': 0.0, 'vary': False}}
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config)
        output_file_name = os.path.join(f'{test_dir}', 'test_wavelength_step_gauss.nxs')
        assert os.path.isfile(output_file_name), f'Expected output file {output_file_name} does not exist.'

        # verify_reduced_data
        gold_file = os.path.join(gold_file_dir, 'expected_wavelength_step_gauss.nxs')
        exp_file = output_file_name
        verify_reduction(exp_file, gold_file, 'gauss')


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
        assert gold_ws.readY(0).shape == test_ws.readY(0),\
            f'Number of wavelength: {gold_ws.readY(0).shape} != {test_ws.readY(0)}'
        assert gold_ws.readX(0).shape == test_ws.readX(0),\
            f'Histogram or point data: {gold_ws.readX(0).shape} != {test_ws.readX(0)}'
        assert np.testing.assert_allclose(gold_ws.extractX(), test_ws.extractX())
        assert np.testing.assert_allclose(gold_ws.extractY(), test_ws.extractY())
        assert np.testing.assert_allclose(gold_ws.extractE(), test_ws.extractE())


if __name__ == '__main__':
    pytest.main([__file__])
