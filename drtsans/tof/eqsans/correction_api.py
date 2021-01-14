# This module contains workflow algorithms and methods correct intensity and error of
# sample and background data accounting wavelength-dependent incoherent inelastic scattering.
# The workflow algorithms will be directly called by eqsans.api.
# from drtsans.tof.eqsans.reduction_api import process_sample_configuration
from drtsans.tof.eqsans.elastic_reference_normalization import (determine_reference_wavelength_q1d_mesh,
                                                                reshape_q_wavelength_matrix,
                                                                determine_common_mod_q_range_mesh,
                                                                calculate_scale_factor_mesh_grid,
                                                                normalize_intensity_q1d,
                                                                build_i_of_q1d)
from collections import namedtuple
"""

Workflow to correct intensities and errors accounting wavelength dependent
incoherent inelastic scattering

Case 1:  If the elastic reference run is NOT specified, the
 - sample
 - bkgd
will be
 * not be corrected in step 2 described in story #689.
 * corrected in step 3 and step 4 described in story #689.

Case 2: If an elastic reference run is specified.
        Step 2 in story #689 will be executed.

  Step 2:
    a) parse additional configuration

        elastic_ref = reduction_input["elasticReference"]["runNumber"]
        elastic_ref_trans = reduction_input["elasticReferenceTrans"]["runNumber"]
        elastic_ref_trans_value = reduction_input["elasticReferenceTrans"]["value"]
        elastic_ref_bkgd = None  # in story 689
        elastic_ref_bkgd_trans = None
        elastic_ref_bkgd_trans_value = None

    b) load elastic reference run, reference transition run, background run and background transition run
       This shall be implemented in eqsans.api.load_all_files()

    c) calculate K and delta K

    d) before bin_all() called in method reduce_single_configuration
       > # Bin 1D and 2D
       > iq2d_main_out, iq1d_main_out = bin_all(iq2d_main_in_fr[wl_frame], iq1d_main_in_fr[wl_frame],
       > ... ...

       1. normalize sample and bkgd run
       2. execute step 3 and step 4
"""


# TODO - when python is upgraded to 3.7+, this class shall be wrapped as dataclass
class CorrectionConfiguration:
    """
    A data class/structure to hold the parameters configured to do incoherence/inelastic
    scattering correction
    """
    def __init__(self, do_correction=False, select_min_incoherence=False):

        self._do_correction = do_correction
        self._select_min_incoherence = select_min_incoherence
        self._elastic_ref_run_setup = None

    @property
    def do_correction(self):
        return self._do_correction

    @property
    def select_min_incoherence(self):
        return self._select_min_incoherence

    @select_min_incoherence.setter
    def select_min_incoherence(self, flag):
        self._select_min_incoherence = flag

    @property
    def elastic_reference_run(self):
        return self._elastic_ref_run_setup

    def set_elastic_reference_run(self, reference_run_setup):
        """Set elastic reference run reduction setup

        Parameters
        ----------
        reference_run_setup: ElasticReferenceRunSetup
            reduction setup

        Returns
        -------

        """
        self._elastic_ref_run_setup = reference_run_setup


# TODO - when python is upgraded to 3.7+, this class shall be wrapped as dataclass
class ElasticReferenceRunSetup:
    """
    A data class/structure to hold the reference run
    """
    def __init__(self, ref_run_number, trans_run_number=None, trans_value=None):
        self.ref_run_number = ref_run_number
        self.trans_run_number = trans_run_number
        self.trans_value = trans_value


def parse_correction_config(reduction_config):
    """Parse correction configuration from reduction configuration

    Parameters
    ----------
    reduction_config: ~dict
        reduction configuration from JSON

    Returns
    -------
    CorrectionConfiguration
        incoherence/inelastic scattering correction configuration
    """
    run_config = reduction_config['configuration']
    do_correction = run_config.get('fitInelasticIncoh', False)
    select_min_incoherence = run_config.get('selectMinIncoh', False)
    _config = CorrectionConfiguration(do_correction, select_min_incoherence)
    elastic_ref = run_config.get('elasticReference')
    if elastic_ref is not None:
        elastic_ref = ElasticReferenceRunSetup(elastic_ref)
        _config.set_elastic_reference_run(elastic_ref)
    return _config


def process_bin_workspace(raw_ws, transmission, ref_sample_thickness, binning_setup):
    """Process rar workspace and do the binning with keeping wavelength terms

    Parameters
    ----------
    raw_ws
    transmission
    ref_sample_thickness
    binning_setup

    Returns
    -------
    ~tuple
        list of IQmod, list of IQazimuthal

    """

    ref_trans_ws, ref_trans_value = transmission
    assert ref_trans_ws
    assert ref_trans_value
    assert raw_ws
    assert ref_sample_thickness
    assert binning_setup

    return ['IofQ1D'], ['IofQ2D']

NormFactor = namedtuple('NormFactor', 'k k_error p s')


def calculate_elastic_scattering_factor(ref_ws, ref_trans_ws, ref_trans_value, ref_sample_thickness,
                                        binning_setup):
    """Normalize runs (sample and background) for elastic scattering

    for elastic reference, I assume that it will use the same as the sample run.
    Then as the elastic reference run is reduced, K and delta K are calculated.
    These runs will be normalized by K and delta K value
    - sample
    - bkgd


    Parameters
    ----------
    ref_ws
    ref_trans_ws
    ref_trans_value
    ref_sample_thickness
    binning_setup

    Returns
    -------
    ~dict
        blabla
    """
    # Process elastic reference run, reference transmission run,
    # TODO - reference background run, reference background transmission run
    ref_iq1d_frames, ref_iq2d_frames = process_bin_workspace(ref_ws, (ref_trans_ws, ref_trans_value),
                                                             ref_sample_thickness, binning_setup)

    # Sanity check
    assert len(ref_iq1d_frames) <= 3, f'Number of frames {len(ref_iq1d_frames)} is not reasonable.'

    # Output
    elastic_norm_factor_dict = dict()

    # Calculate scaling vectors for each
    for i_frame in len(ref_iq1d_frames):
        # reshape
        ref_wl_vec, ref_q_vec, ref_i_array, ref_error_array, ref_dq_array = reshape_q_wavelength_matrix(
            ref_iq1d_frames[i_frame])

        # Calculate Qmin and Qmax
        qmin_index, qmax_index = determine_common_mod_q_range_mesh(ref_q_vec, ref_i_array)

        # Calculate reference
        ref_wl_ie = determine_reference_wavelength_q1d_mesh(ref_wl_vec, ref_q_vec, ref_i_array, ref_error_array,
                                                            qmin_index, qmax_index)

        # Calculate scale factor
        k_vec, k_error_vec, p_vec, s_vec = calculate_scale_factor_mesh_grid(ref_wl_vec, ref_i_array, ref_error_array,
                                                                            ref_wl_ie, qmin_index, qmax_index)

        norm_factor = NormFactor(k_vec, k_error_vec, p_vec, s_vec)
        # Form output
        elastic_norm_factor_dict[i_frame] = norm_factor

        # export K and K error (task #725)
        # do_export_k(k_vec, k_error_vec, k_filename)

    return elastic_norm_factor_dict


def normalize_ws_with_elastic_scattering(i_q1d_frames, i_q2d_frames, norm_dict):

    # KEY: a data structure to contain IQmod(s) of (1) all samples and background (2) all frames
    if i_q2d_frames:
        print(f'Not implemented!')

    # Normalize sample and background
    # normalize 1D
    num_frames = len(i_q1d_frames)
    norm_iq1d = list()
    for i_frame in range(num_frames):
        # convert
        wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_q1d_frames[i_frame])
        qmin_index = qmax_index = -1
        if qmin_index + qmax_index < 0:
            raise RuntimeError('calcualte qmin and qmax!')
        # reference wavelength
        data_ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                                 qmin_index, qmax_index)
        #
        k_vec = norm_dict[i_frame].k
        p_vec = norm_dict[i_frame].p
        s_vec = norm_dict[i_frame].s

        # normalize
        normalized = normalize_intensity_q1d(wl_vec, q_vec, i_array, error_array,
                                             data_ref_wl_ie, k_vec, p_vec, s_vec,
                                             qmin_index, qmax_index)

        # Convert normalized intensities and errors to IModQ
        normalized_i_of_q = build_i_of_q1d(wl_vec, q_vec, normalized[0], normalized[1], dq_array)

        # set
        norm_iq1d.append(normalized_i_of_q)

    return norm_iq1d


# def process_elastic_runs():
#     """
#     - dark current
#     - mask
#     - sensitivities
#     - empty beam
#     - empty beam radius
#     Returns
#     -------
#
#     """
#     processed_elastic_ref = process_sample_configuration(elastic_reference_run,
#                                                          sample_trans_ws=elastic_reference_ws,
#                                                          sample_trans_value=elastic_reference_ws_value,
#                                                          bkg_ws_raw=None,
#                                                          bkg_trans_ws=None,  # bkgd_trans_ws,
#                                                          bkg_trans_value=None,  # bkg_trans_value,
#                                                          theta_deppendent_transmission=True / False,
#                                                          dark_current=loaded_ws.dark_current,
#                                                          flux_method=flux_method,
#                                                          flux=flux,
#                                                          mask_ws=loaded_ws.mask,
#                                                          mask_panel=mask_panel,
#                                                          solid_angle=solid_angle,
#                                                          sensitivity_workspace=loaded_ws.sensitivity,
#                                                          output_workspace=f'processed_data_main',
#                                                          output_suffix=output_suffix,
#                                                          thickness=thickness,
#                                                          absolute_scale_method=absolute_scale.method,
#                                                          # absolute_scale_method,
#                                                          absolute_scale=absolute_scale.value,  # absolute_scale,
#                                                          empty_beam_ws=empty_trans_ws,
#                                                          beam_radius=beam_radius,
#                                                          keep_processed_workspaces=False)
#
#     return i_of_q


def correct_acc_incoherence_scattering(iq1d_frames, iq2d_frames, correction_setup):
    assert iq1d_frames
    assert iq2d_frames
    assert correction_setup

    return iq1d_frames, iq2d_frames
