# This module contains workflow algorithms and methods correct intensity and error of
# sample and background data accounting wavelength-dependent incoherent inelastic scattering.
# The workflow algorithms will be directly called by eqsans.api.
from drtsans.tof.eqsans.reduction_api import process_sample_configuration
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


from drtsans.tof.eqsans.reduction_api import process_sample_configuration

def correct_main_workflow():

    # process elastic reference run
    process_sample_configuration(elastic_run, elastic_trans, elastic_trans_value)




def calculate_elastic_scattering_factor(runs, ref_thickness):
    # normalize runs (sample and background) for elastic scattering
    """
    for elastic reference, I assume that it will use the same as the sample run.

    Then as the elastic reference run is reduced, K and delta K are calculated.

    These runs will be normalized by K and delta K value
    - sample
    - bkgd

    Parameters
    ----------
    runs: ~list
        sample and background run (workspace)
    ref_thickness: float
        elastic reference sample thickness in unit of mm

    Returns
    -------
    ~list

    """
    # process elastic run
    ref_i_of_q_frames = process_elastic_runs()

    # calculate K and K error
    from drtsans.tof.eqsans.elastic_reference_normalization import (determine_reference_wavelength_q1d_mesh,
                                                                    reshape_q_wavelength_matrix,
                                                                    determine_common_mod_q_range_mesh,
                                                                    calculate_scale_factor_mesh_grid,
                                                                    normalize_intensity_q1d,
                                                                    build_i_of_q1d)

    # sanity check
    assert len(ref_i_of_q_frames) <= 3, 'Shape in correct!'
    for i_frame in len(ref_i_of_q_frames):
        # reshape
        ref_wl_vec, ref_q_vec, ref_i_array, ref_error_array, ref_dq_array = reshape_q_wavelength_matrix(
            ref_i_of_q_frames[i_frame])

        # Calculate Qmin and Qmax
        qmin_index, qmax_index = determine_common_mod_q_range_mesh(ref_q_vec, ref_i_array)

        # Calculate reference
        ref_wl_ie = determine_reference_wavelength_q1d_mesh(ref_wl_vec, ref_q_vec, ref_i_array, ref_error_array,
                                                            qmin_index, qmax_index)

        # Calculate scale factor
        k_vec, k_error_vec, p_vec, s_vec = calculate_scale_factor_mesh_grid(ref_wl_vec, ref_i_array, ref_error_array,
                                                                            ref_wl_ie, qmin_index, qmax_index)

        #
        k_vec_set[i_frame] = k_vec
        k_error_set[i_frame] = k_error_vec

        # export K and K error (task #725)
        # do_export_k(k_vec, k_error_vec, k_filename)

    return k_vec_set, k_error_set


def normalize_elastic_scattering():

    # KEY: a data structure to contain IQmod(s) of (1) all samples and background (2) all frames

    # Normalize sample and background
    # normalize 1D
    normalized_set = list()
    for i_of_q in [combined_i_of_q(sample_i_of_q_list, bkgd_i_of_q)]:
        # loop over frames to normalize
        for i_frame in range(num_frames):
            wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_of_q[i_frame])
            data_ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                                     qmin_index, qmax_index)
        normalized = normalize_intensity_q1d(wl_vec, q_vec, i_array, error_array,
                                             data_ref_wl_ie, k_vec, p_vec, s_vec,
                                             qmin_index, qmax_index)

        # Convert normalized intensities and errors to IModQ
        normalized_i_of_q = build_i_of_q1d(wl_vec, q_vec, normalized[0], normalized[1], dq_array)
        normalized_set.append(normalized_i_of_q)

    return normalized_set[0], normalized_set[1]


def process_elastic_runs():
    """
    - dark current
    - mask
    - sensitivities
    - empty beam
    - empty beam radius
    Returns
    -------

    """
    processed_elastic_ref = process_sample_configuration(elastic_reference_run,
                                                         sample_trans_ws=elastic_reference_ws,
                                                         sample_trans_value=elastic_reference_ws_value,
                                                         bkg_ws_raw=None,
                                                         bkg_trans_ws=None,  # bkgd_trans_ws,
                                                         bkg_trans_value=None,  # bkg_trans_value,
                                                         theta_deppendent_transmission=True / False,
                                                         dark_current=loaded_ws.dark_current,
                                                         flux_method=flux_method,
                                                         flux=flux,
                                                         mask_ws=loaded_ws.mask,
                                                         mask_panel=mask_panel,
                                                         solid_angle=solid_angle,
                                                         sensitivity_workspace=loaded_ws.sensitivity,
                                                         output_workspace=f'processed_data_main',
                                                         output_suffix=output_suffix,
                                                         thickness=thickness,
                                                         absolute_scale_method=absolute_scale.method,
                                                         # absolute_scale_method,
                                                         absolute_scale=absolute_scale.value,  # absolute_scale,
                                                         empty_beam_ws=empty_trans_ws,
                                                         beam_radius=beam_radius,
                                                         keep_processed_workspaces=False)

    return i_of_q
