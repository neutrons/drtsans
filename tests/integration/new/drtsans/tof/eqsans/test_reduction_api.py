import pytest
import os
import tempfile
from drtsans.tof.eqsans import reduction_parameters, update_reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch
import numpy as np
from drtsans.dataobjects import save_i_of_q_to_h5, load_iq1d_from_h5, load_iq2d_from_h5
from typing import List, Any, Union
from drtsans.dataobjects import _Testing


# EQSANS reduction
specs_eqsans = {
    'EQSANS_88980': {
        "iptsNumber": 19800,
        "sample": {"runNumber": 88980, "thickness": 0.1, "transmission": {"runNumber": 88980}},
        "background": {"runNumber": 88978, "transmission": {"runNumber": 88974}},
        "beamCenter": {"runNumber": 88973},
        "emptyTransmission": {"runNumber": 88973},
        "configuration": {
            "sampleApertureSize": 30,
            "darkFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5",
            "StandardAbsoluteScale": 0.0208641883,
            "sampleOffset": 0,
        }
    }
}


@pytest.mark.parametrize('run_config, basename',
                         [(specs_eqsans['EQSANS_88980'], 'EQSANS_88980')],
                         ids=['88980'])
def test_regular_setup(run_config, basename, tmpdir, reference_dir):
    """Same reduction from Shaman test

    Returns
    -------

    """
    common_config = {
        "configuration": {
            "maskFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs",
            "useDefaultMask": True,
            "normalization": "Total charge",
            "fluxMonitorRatioFile": "/SNS/EQSANS/IPTS-24769/shared/EQSANS_110943.out",
            "beamFluxFileName": "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample",
            "absoluteScaleMethod": "standard",
            "detectorOffset": 0,
            "mmRadiusForTransmission": 25,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 120,
            "AnnularAngleBin": 5,
            "wavelengthStepType": "constant Delta lambda",
            "wavelengthStep": 0.1,
        }
    }
    input_config = reduction_parameters(common_config, 'EQSANS', validate=False)  # defaults and common options
    input_config = update_reduction_parameters(input_config, run_config, validate=False)
    output_dir = str(tmpdir)
    amendments = {
        'outputFileName': basename,
        'configuration': {'outputDir': output_dir}
    }
    input_config = update_reduction_parameters(input_config, amendments, validate=True)  # final changes and validation

    # Load and reduce
    loaded = load_all_files(input_config)
    reduction_output = reduce_single_configuration(loaded, input_config)

    # Load data and compare
    gold_dir = reference_dir.new.eqsans
    for index in range(2):
        # 1D
        iq1d_h5_name = os.path.join(gold_dir, f'gold_iq1d_{index}_0.h5')
        gold_iq1d = load_iq1d_from_h5(iq1d_h5_name)
        _Testing.assert_allclose(reduction_output[index].I1D_main[0], gold_iq1d)

        # 2D
        iq2d_h5_name = os.path.join(gold_dir, f'gold_iq2d_{index}.h5')
        gold_iq2d = load_iq2d_from_h5(iq2d_h5_name)
        _Testing.assert_allclose(reduction_output[index].I2D_main, gold_iq2d)


def export_reduction_output(reduction_output: List[Any], output_dir: Union[None, str] = None, prefix: str = ''):
    """Export the reduced I(Q) and I(Qx, Qy) to  hdf5 files
    """
    # FIXME TODO - need to document to eqsans.api
    # output of reduce_single_configuration: list of list, containing
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
        reduction_output = reduce_single_configuration(loaded, input_config)
        output_file_name = os.path.join(test_dir, 'test_wavelength_step_reg.nxs')
        assert os.path.isfile(output_file_name), f'Expected output file {output_file_name} does not exists'

        # verify_reduced_data
        gold_file = os.path.join(gold_file_dir, 'expected_wavelength_step_reg.nxs')
        exp_file = output_file_name
        verify_reduction(exp_file, gold_file, 'reg')
        # verify binned IQmod and IQazimuthal
        gold_dir = reference_dir.new.eqsans
        # 1D
        iq1d_h5_name = os.path.join(gold_dir, f'gold_iq1d_wave_0_0.h5')
        gold_iq1d = load_iq1d_from_h5(iq1d_h5_name)
        _Testing.assert_allclose(reduction_output[0].I1D_main[0], gold_iq1d)

        # 2D
        iq2d_h5_name = os.path.join(gold_dir, f'gold_iq2d_wave_0.h5')
        gold_iq2d = load_iq2d_from_h5(iq2d_h5_name)
        _Testing.assert_allclose(reduction_output[0].I2D_main, gold_iq2d)

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


if __name__ == '__main__':
    pytest.main([__file__])
