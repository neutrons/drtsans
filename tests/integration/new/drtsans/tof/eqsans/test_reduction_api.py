import pytest
import os
import tempfile
from drtsans.tof.eqsans import reduction_parameters, update_reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from mantid.simpleapi import LoadNexusProcessed, CheckWorkspacesMatch
import numpy as np
from drtsans.dataobjects import save_i_of_q_to_h5, load_iq1d_from_h5, load_iq2d_from_h5
from typing import List, Any, Union, Tuple
from drtsans.dataobjects import _Testing
from matplotlib import pyplot as plt
from drtsans.dataobjects import IQmod


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
def test_correction_workflow(run_config, basename, tmpdir, reference_dir):
    """Same reduction from Shaman test but using the workflow that is designed to work with inelastic correction

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
            "useErrorWeighting": True,
        }
    }
    input_config = reduction_parameters(common_config, 'EQSANS', validate=False)  # defaults and common options
    input_config = update_reduction_parameters(input_config, run_config, validate=False)
    output_dir = str(tmpdir)
    amendments = {
        'outputFileName': f'{basename}_corr',
        'configuration': {'outputDir': output_dir}
    }
    input_config = update_reduction_parameters(input_config, amendments, validate=True)  # final changes and validation

    # expected output Nexus file
    reduced_data_nexus = os.path.join(output_dir, f'{basename}_corr.nxs')
    # clean output directory
    if os.path.exists(reduced_data_nexus):
        os.remove(reduced_data_nexus)

    # Test: use correction workflow and remove background
    use_correction_workflow = True  # True: new workflow to test; False: old workflow
    keep_background = False  # default to be False

    # Load and reduce
    loaded = load_all_files(input_config)
    reduction_output = reduce_single_configuration(loaded, input_config,
                                                   use_correction_workflow=use_correction_workflow,
                                                   ignore_background=keep_background)

    gold_dir = reference_dir.new.eqsans
    # Verify existence of reduced
    reduced_data_nexus = os.path.join(output_dir, f'{basename}_corr.nxs')
    assert os.path.exists(reduced_data_nexus), f'Expected {reduced_data_nexus} does not exist'
    # Verify reduced workspace (previous) FIXME - remove later
    gold_ws_nexus = os.path.join(reference_dir.new.eqsans, 'EQSANS_88980_reduced.nxs')
    verify_reduction(test_file=reduced_data_nexus,  gold_file=gold_ws_nexus, ws_prefix='no_wl', ignore_error=True)
    # verify with gold data
    gold_ws_nexus = os.path.join(gold_dir, 'test_integration_api/EQSANS_88980_reduced_wb_m6.nxs')
    verify_reduction(test_file=reduced_data_nexus,  gold_file=gold_ws_nexus, ws_prefix='no_wl', ignore_error=False)

    # Verify binned I(Q)
    for index in range(2):
        # Load expected I(Q) and I(Qx, Qy)
        # Fixme : these 'expected' intensity and error are not verified but only Q range
        gold_iq1d_h5 = os.path.join(gold_dir, f'88980_frame1_weighted_old_removebkgd_{index}.h5')
        assert os.path.exists(gold_iq1d_h5)
        gold_iq1d = load_iq1d_from_h5(gold_iq1d_h5)

        gold_iq2d_h5 = os.path.join(gold_dir, f'gold_iq2d_{index}.h5')
        assert os.path.exists(gold_iq2d_h5)
        gold_iq2d = load_iq2d_from_h5(gold_iq2d_h5)
        print(f'Verifying intensity frame {index} from {gold_iq1d_h5} and {gold_iq2d_h5}')

        # Verify Q bins: 1D only, 2D skip
        np.testing.assert_allclose(gold_iq1d.mod_q, reduction_output[index].I1D_main[0].mod_q)
        # 2D
        assert gold_iq2d

    """ This is information about how gold data will be generated for the next step: binning
    on different reduction workflow
    # Save I(Q) to h5 file
    flag1 = 'old' if not use_correction_workflow else 'new'
    flag2 = 'keepbkgd' if keep_background else 'removebkgd'
    for index in range(2):
        save_i_of_q_to_h5(reduction_output[index].I1D_main[0],
                          os.path.join(cwd, f'88980_frame1_weighted_{flag1}_{flag2}_{index}.h5'))
    gold_iq1d_h5 = os.path.join(gold_dir, f'88980_frame1_weighted_old_removebkgd_{index}.h5')
    gold_iq2d_h5 = os.path.join(gold_dir, f'gold_iq2d_{index}.h5')
    """


@pytest.mark.parametrize('run_config, basename',
                         [(specs_eqsans['EQSANS_88980'], 'EQSANS_88980')],
                         ids=['88980'])
def test_regular_setup(run_config, basename, tmpdir, reference_dir):
    """Same reduction from Shaman test

    Returns
    -------

    """
    # set flag to use weighted binning
    weighted_binning = False

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
            "useErrorWeighting": weighted_binning,
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

    # expected output Nexus file
    reduced_data_nexus = os.path.join(output_dir, f'{basename}.nxs')
    # remove files
    if os.path.exists(reduced_data_nexus):
        os.remove(reduced_data_nexus)

    # Load and reduce
    loaded = load_all_files(input_config)
    reduction_output = reduce_single_configuration(loaded, input_config, use_correction_workflow=False)

    # Load data and compare
    gold_dir = reference_dir.new.eqsans
    for index in range(2):
        # 1D
        iq1d_h5_name = os.path.join(gold_dir, f'test_integration_api/88980_iq1d_{index}_0_m6.h5')
        gold_iq1d = load_iq1d_from_h5(iq1d_h5_name)
        export_iq_comparison([('Test Result', reduction_output[index].I1D_main[0], 'red'),
                              ('Gold Result', gold_iq1d, 'green')],
                             f'/tmp/regular_setup_comparison_{index}.png')
        _Testing.assert_allclose(reduction_output[index].I1D_main[0], gold_iq1d)

        # 2D
        iq2d_h5_name = os.path.join(gold_dir, f'test_integration_api/88980_iq2d_{index}_m6.h5')
        gold_iq2d = load_iq2d_from_h5(iq2d_h5_name)
        _Testing.assert_allclose(reduction_output[index].I2D_main, gold_iq2d)

    # Check reduced workspace
    assert os.path.exists(reduced_data_nexus), f'Expected {reduced_data_nexus} does not exist'
    # verify with gold data and clean
    gold_file = os.path.join(reference_dir.new.eqsans, 'test_integration_api/EQSANS_88980_reduced_m6.nxs')
    verify_reduction(test_file=reduced_data_nexus,  gold_file=gold_file, ws_prefix='no_wl')
    # clean up
    os.remove(reduced_data_nexus)


def export_iq_comparison(iq1d_tuple_list: List[Tuple[str, IQmod, str]], png_name: str):
    """Export a list of IQmod to plot
    """

    plt.figure(figsize=(18, 9))
    for iq1d_tuple in iq1d_tuple_list:
        label, iq1d, color = iq1d_tuple
        plt.plot(iq1d.mod_q, iq1d.intensity, color=color, label=label)

    # legend
    plt.legend()

    # save
    plt.savefig(png_name)
    # close
    plt.close()

    plt.figure(figsize=(18, 12))
    plt.yscale('log')

    # plot error bar to compare
    for iq1d_tuple in iq1d_tuple_list:
        label, iq1d, color = iq1d_tuple
        plt.plot(iq1d.mod_q, iq1d.error, color=color, label=label, marker='.', linestyle='None')

    # legend
    plt.legend()

    # save
    plt.savefig(f'{png_name.split(".")[0]}_error_bar.png')
    # close
    plt.close()


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

    # Specify gold dir
    gold_dir = reference_dir.new.eqsans

    # Test 1 with regular setup
    with tempfile.TemporaryDirectory() as test_dir:
        # continue to configure
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_reg'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        # reduce
        loaded = load_all_files(input_config)
        reduction_output = reduce_single_configuration(loaded, input_config)

        # verify output file existence
        output_file_name = os.path.join(test_dir, 'test_wavelength_step_reg.nxs')
        assert os.path.isfile(output_file_name), f'Expected output file {output_file_name} does not exists'
        # verify reduced worksapce
        gold_file = os.path.join(gold_dir, 'test_integration_api/EQSANS_88980_wl_reduced_reg_m6.nxs')
        verify_reduction(output_file_name, gold_file, 'reg', ignore_error=False)
        # verify binned reduced I(Q)
        gold_iq1d_h5 = os.path.join(gold_dir, 'test_integration_api/88980_iq1d_wl_0_0_m6.h5')
        gold_iq1d = load_iq1d_from_h5(gold_iq1d_h5)
        _Testing.assert_allclose(reduction_output[0].I1D_main[0], gold_iq1d)
        # verify binned reduced I(Qx, Qy)
        # TODO skip as no knowing what the user's requirement with wavelength kept
        # iq2d_h5_name = os.path.join(gold_dir, f'gold_iq2d_wave_0.h5')
        # gold_iq2d = load_iq2d_from_h5(iq2d_h5_name)
        # test_iq2d = reduction_output[0].I2D_main
        # _Testing.assert_allclose(reduction_output[0].I2D_main, gold_iq2d)

    # Test 2 with c.o.m beam center
    with tempfile.TemporaryDirectory() as test_dir:
        # continue to configure
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_com'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        input_config['beamCenter']['method'] = 'center_of_mass'
        input_config['beamCenter']['com_centering_options'] = {'CenterX': 0., 'CenterY': 0., 'Tolerance': 0.00125}
        # reduce
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config)
        # verify output file existence
        output_file_name = os.path.join(f'{test_dir}', 'test_wavelength_step_com.nxs')
        assert os.path.isfile(output_file_name),  f'Expected output file {output_file_name} does not exists'
        # verify reduced worksapce
        gold_file = os.path.join(gold_dir, 'test_integration_api/EQSANS_88980_wl_reduced_com_m6.nxs')
        verify_reduction(output_file_name, gold_file, 'com', ignore_error=False)

    # Test 3 with gaussian beam center
    with tempfile.TemporaryDirectory() as test_dir:
        # continue to configure
        configuration['configuration']['outputDir'] = test_dir
        configuration['outputFileName'] = 'test_wavelength_step_gauss'
        configuration['dataDirectories'] = test_dir
        # validate and clean configuration
        input_config = reduction_parameters(configuration)
        input_config['beamCenter']['method'] = 'gaussian'
        input_config['beamCenter']['gaussian_centering_options'] = {'theta': {'value': 0.0, 'vary': False}}
        # reduce
        loaded = load_all_files(input_config)
        reduce_single_configuration(loaded, input_config)

        # verify output file existence
        output_file_name = os.path.join(f'{test_dir}', 'test_wavelength_step_gauss.nxs')
        assert os.path.isfile(output_file_name), f'Expected output file {output_file_name} does not exist.'
        # verify_reduced_data
        # Difference from mantid5 result
        # E   AssertionError:
        # E   Not equal to tolerance rtol=1e-07, atol=0   Y is not same
        # E   Mismatched elements: 1 / 442368 (0.000226%)
        # E   Max absolute difference: 2.96006469e-09  Max relative difference: 1.7555871e-07
        gold_file = os.path.join(gold_dir, 'test_integration_api/EQSANS_88980_wl_reduced_gauss_m6.nxs')
        verify_reduction(output_file_name, gold_file, 'gauss', ignore_error=False, y_rel_tol=3.E-7)


def verify_reduction(test_file, gold_file, ws_prefix, ignore_error=False, y_rel_tol=None):
    """Verify reduced result by verified expected result

    Parameters
    ----------
    test_file: str
        NeXus file from test to verify
    gold_file: str
        NeXus file containing the expected reduced result to verify against
    ws_prefix: str
        prefix for Mantid workspace that the
    ignore_error: bool
        flag to ignore the checking on intensity error
    y_rel_tol: float, None
        allowed maximum tolerance on Y

    """
    assert os.path.exists(gold_file), f'Gold file {gold_file} cannot be found'
    assert os.path.exists(test_file), f'Test file {test_file} cannot be found'

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
        np.testing.assert_allclose(gold_ws.extractX(), test_ws.extractX(), err_msg='X is not same')
        if y_rel_tol is not None:
            y_dict = {'rtol': y_rel_tol}
        else:
            y_dict = dict()
        np.testing.assert_allclose(gold_ws.extractY(), test_ws.extractY(), err_msg='Y is not same', **y_dict)
        if not ignore_error:
            np.testing.assert_allclose(gold_ws.extractE(), test_ws.extractE(), err_msg='E is not same')


if __name__ == '__main__':
    pytest.main([__file__])
