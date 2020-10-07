import pytest
import tempfile
import h5py
import numpy as np
import os
from drtsans.mono.spice_data import map_to_nexus


def get_iq1d(log_file_name):
    """Get I(Q) from output SANS log file

    Parameters
    ----------
    log_file_name: str
        log file name

    Returns
    -------
    tuple
        numpy 1D array for Q, numpy 1D array for intensity

    """
    # Open file and entry
    log_h5 = h5py.File(log_file_name, 'r')

    if '_slice_1' in log_h5:
        data_entry = log_h5['_slice_1']['main']
    else:
        data_entry = log_h5['main']

    # Get data
    iq1d_entry = data_entry['I(Q)']

    # Get data with a copy
    vec_q = np.copy(iq1d_entry['Q'].value)
    vec_i = np.copy(iq1d_entry['I'].value)

    # close file
    log_h5.close()

    return vec_q, vec_i


def compare_reduced_iq(test_log_file, gold_log_file, title, prefix):
    """Compare I(Q) from reduced file and gold file

    Parameters
    ----------
    test_log_file
    gold_log_file
    title: str
        title of output figure
    prefix: str
        prefix of output file

    Returns
    -------

    """
    # Plot main
    test_q_vec, test_intensity_vec = get_iq1d(test_log_file)
    gold_q_vec, gold_intensity_vec = get_iq1d(gold_log_file)

    # Verify result
    try:
        np.testing.assert_allclose(test_q_vec, test_q_vec, atol=1E-4)
        np.testing.assert_allclose(test_intensity_vec, gold_intensity_vec, atol=1E-7)
    except AssertionError as assert_err:
        from matplotlib import pyplot as plt
        plt.cla()
        plt.plot(test_q_vec, test_intensity_vec, color='red', label='Corrected')
        plt.plot(gold_q_vec, gold_intensity_vec, color='black', label='Before being corrected')
        plt.legend()
        plt.title(title)
        plt.yscale('log')
        out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '.png'
        plt.savefig(out_name)

        raise assert_err


def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix):

    unmatched_errors = ''

    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        title_i = '{}: {}'.format(sample_name, title)
        try:
            compare_reduced_iq(output_log_file, gold_log_file, title_i, prefix)
        except AssertionError as unmatched_error:
            unmatched_errors = 'Testing output {} is different from gold result {}:\n{}' \
                               ''.format(output_log_file, gold_log_file, unmatched_error)
    # END-FOR

    # raise error for all
    if unmatched_errors != '':
        raise AssertionError(unmatched_errors)


def test_reduction_spice(reference_dir, cleanfile):
    """
    Test reduction from SPICE-converted Nexus file

    """
    CG2 = 'CG2'
    nexus_dir = os.path.join(reference_dir.new.gpsans, 'Exp280')

    # Set output (temp) directory
    output_directory = tempfile.mkdtemp(prefix='spice_reduction')
    # TODO cleanfile(output_directory)

    # USER INPUT
    ipts_number = 828
    exp_number = 280

    # single sample
    samples = [(35, 1)]
    samples_trans = [(27, 1)]
    sample_thick = ['0.1']
    sample_names = ['Porasil_B']
    bkgd = [(34, 1)]
    bkgd_trans = [(26, 1)]

    empty_trans = [(28, 1)]
    beam_center = [(20, 1)]

    # q range to use to clean 1D curve of each configuration
    q_range = [None, None]

    # STAFF INPUT
    use_mask_file = True
    mask_file_name = os.path.join(reference_dir.new.gpsans, 'calibrations/mask_pixel_map.nxs')
    use_dark_file = False
    dark_file_name = ""
    block_beam = (9, 1)
    use_mask_back_tubes = False
    wavelength = None
    wavelength_spread = None
    wedge_min_angles = None
    wedge_max_angles = None

    # sensitivity_file = '/HFIR/CG2/shared/drt_sensitivity/sens_CG2_spice_bar.nxs'
    sensitivity_file = os.path.join(reference_dir.new.gpsans, 'calibrations/sens_CG2_spice_bar.nxs')
    use_log_2d_binning = False
    use_log_1d = True
    common_configuration = {
        "iptsNumber": ipts_number,
        "emptyTransmission": {"runNumber": map_to_nexus(CG2, ipts_number, exp_number, empty_trans, nexus_dir)},
        "beamCenter": {"runNumber": map_to_nexus(CG2, ipts_number, exp_number, beam_center, nexus_dir)},
        "configuration": {
            "outputDir": output_directory,
            "darkFileName": dark_file_name,
            "sensitivityFileName": sensitivity_file,
            "DBScalingBeamRadius": 40,
            "sampleApertureSize": 8,
            "mmRadiusForTransmission": 40,
            "absoluteScaleMethod": "direct_beam",
            "numQxQyBins": 256,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numQBins": "",
            "LogQBinsPerDecade": 33,
            "useLogQBinsDecadeCenter": True,
            "useLogQBinsEvenDecade": False,
            "wavelength": wavelength,
            "wavelengthSpread": wavelength_spread,
            "blockedBeamRunNumber": 'Whatever',
            "maskFileName": mask_file_name,
            'WedgeMinAngles': wedge_min_angles,
            'WedgeMaxAngles': wedge_max_angles,
            "AnnularAngleBin": 2.0,
            "Qmin": 0.0028,
            "Qmax": 0.0035,
            "useSubpixels": True,
            "subpixelsX": 5,
            "subpixelsY": 5,
            "useTimeSlice": False,
            "useLogSlice": False,
            "logSliceName": "",
            "logSliceInterval": '',
        }
    }

    # Never touch!  drtsans specific

    # convert SPICE to Nexus
    samples = map_to_nexus(CG2, ipts_number, exp_number, samples, nexus_dir)
    samples_trans = map_to_nexus(CG2, ipts_number, exp_number, samples_trans, nexus_dir)
    bkgd = map_to_nexus(CG2, ipts_number, exp_number, bkgd, nexus_dir)
    bkgd_trans = map_to_nexus(CG2, ipts_number, exp_number, bkgd_trans, nexus_dir)

    import warnings
    warnings.filterwarnings('ignore')
    import json
    # jupyter only:from pprint import pprint as pretty_print
    from drtsans.mono.gpsans import (load_all_files, reduce_single_configuration, plot_reduction_output,
                                     reduction_parameters, update_reduction_parameters)
    from matplotlib.colors import LogNorm

    if not use_mask_file:
        mask_file_name = ""
    if not use_dark_file:
        dark_file_name = ""

    if use_log_2d_binning:
        log_flag = {"norm": LogNorm()}
    else:
        log_flag = {'vmin': 0, 'vmax': 100}

    # Add on the other reduction parameters with their default values (most will be empty)
    common_configuration_full = reduction_parameters(common_configuration, 'GPSANS', validate=False)

    # Create output directory
    output_dir = common_configuration_full['configuration']['outputDir']
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    for i in range(len(samples)):
        # Settings particular to each reduction session
        run_data = {
            'sample': {
                'runNumber': samples[i],
                'thickness': sample_thick[i],
                'transmission': {'runNumber': samples_trans[i]}
            },
            'background': {
                'runNumber': bkgd[i],
                'transmission': {'runNumber': bkgd_trans[i]}
            },
            'outputFileName': sample_names[i],
            'configuration': {
                "Qmin": q_range[0],
                "Qmax": q_range[1],
                "useMaskBackTubes": use_mask_back_tubes,
                "blockedBeamRunNumber": map_to_nexus(CG2, ipts_number, exp_number, [block_beam], nexus_dir)[0],
                "maskFileName": mask_file_name,
                "darkFileName": dark_file_name,
            }
        }

        # Update our common settings with the particulars of the current reduction
        reduction_input = update_reduction_parameters(common_configuration_full, run_data, validate=True)

        # Begin reduction. Be sure to validate the parameters before.
        loaded = load_all_files(reduction_input, path=f'/HFIR/CG2/IPTS-{ipts_number}/shared/Exp{exp_number}')
        out = reduce_single_configuration(loaded, reduction_input)
        plot_reduction_output(out, reduction_input, loglog=use_log_1d, imshow_kwargs=log_flag)

        # Save the reduction parameters of each reduction session to a JSON file
        output_dir = reduction_input['configuration']['outputDir']
        output_json_file = os.path.join(output_dir, f'{sample_names[i]}.json')  # full path to the JSON file
        with open(output_json_file, 'w') as file_handle:
            json.dump(reduction_input, file_handle, indent=2)

        # verify
        expected_data_dir = os.path.join(reference_dir.new.gpsans, 'reduced_exp280')
        verify_reduction_results(sample_names, output_dir, expected_data_dir, 'SPICE reduction', prefix=None)


if __name__ == '__main__':
    pytest.main([__file__])
