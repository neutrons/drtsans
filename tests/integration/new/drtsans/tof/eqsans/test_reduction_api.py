import pytest
import os
import tempfile
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason="Required test data not available")
def test_wavelength_step():

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
        assert os.path.isfile(output_file_name), f'File {output_file_name} exists = {os.path.exists(output_file_name)}'

        # verify_reduced_data
        gold_file = 'expected_wavelength_step_reg.nxs'
        exp_file = output_file_name
        assert os.path.exists(gold_file)
        gold_ws = LoadNexusProcessed(Filename=gold_file, OutputWorkspace='gold_reg')
        test_ws = LoadNexusProcessed(Filename=exp_file, OutputWorkspace='test_reg')
        r = CheckWorkspacesMatch(Workspace1=gold_ws, Workspace2=test_ws)
        assert r == 'Success!'

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
        assert os.path.isfile(os.path.join(f'{test_dir}', 'test_wavelength_step_com.nxs'))
        # assert os.path.isfile(f'{test_dir}/test_wavelength_step.nxs')

        # verify_reduced_data
        gold_file = 'expected_wavelength_step_com.nxs'
        exp_file = output_file_name
        assert os.path.exists(gold_file)
        gold_ws = LoadNexusProcessed(Filename=gold_file, OutputWorkspace='gold_com')
        test_ws = LoadNexusProcessed(Filename=exp_file, OutputWorkspace='test_com')
        r = CheckWorkspacesMatch(Workspace1=gold_ws, Workspace2=test_ws)
        assert r == 'Success!'

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
        assert os.path.isfile(output_file_name), f'Output file {output_file_name} does not exist.'

        # verify_reduced_data
        gold_file = 'expected_wavelength_step_gauss.nxs'
        exp_file = output_file_name
        assert os.path.exists(gold_file)
        gold_ws = LoadNexusProcessed(Filename=gold_file, OutputWorkspace='gold_guass')
        test_ws = LoadNexusProcessed(Filename=exp_file, OutputWorkspace='test_gauss')
        r = CheckWorkspacesMatch(Workspace1=gold_ws, Workspace2=test_ws)
        assert r == 'Success!'


if __name__ == '__main__':
    pytest.main([__file__])
