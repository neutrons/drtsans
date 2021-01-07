# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
# Step 3: Correct wave-dependent incoherence intensity for I(Q, wavelength)
from drtsans.tof.eqsans.elastic_reference_normalization import determine_common_mod_q_range
from drtsans.tof.eqsans.elastic_reference_normalization import determine_reference_wavelength_q1d
import numpy as np


def correct_incoherence_inelastic_1d(i_of_q):

    # Verify: q-bins are same for all wavelength bins: meshgrid between Q and wavelength
    # TODO: based on the requirement as q-bins must be same for all wavelength bins,
    #       it can be pushed backward to the start of the workflow that Q and wavelength are in mesh grid

    # Determine q_min and q_max  that exist in all wavelength for the fitting (minimization) process.
    min_q, max_q = determine_common_mod_q_range(i_of_q)

    # Calculate reference wavelength for each Q: shortest wavelength bin
    ref_wl_matrix = determine_reference_wavelength_q1d(iqmod)

    # Calculate inelastic incoherent factor: b(lambda)
    b_factor_vec = calculate_inelastic_incoherent_factor(i_of_q, min_q, max_q)

    # Optionally select minimum incoherence ...
    if select_min_incoh:
        # re-select minimum incoherence wavelength
        ref_wl_matrix2 = select_min_incoherence(b_factor_vec, unique_wl)
        # re-alculate inelastic incoherent factor: b(lambda)
        b_factor_vec = calculate_inelastic_incoherent_factor(i_of_q, min_q, max_q, ref_wl_matrix2)
        # check: all b_factor_vec[i_wl] shall be greater than ZERO
        assert b_factor_vec > 0

    # Update data for I(Q, wavelength) and delta I(Q, wavelength)


def calculate_inelastic_incoherent_factor(iqmod, q_min, q_max, ref_wl_matrix):

    # Create a numpy 2D array for slicing/filtering: wavelength, qmod, intensity
    wqi_matrix = np.array([iqmod.wavelength, iqmod.mod_q, iqmod.intensity])
    wqi_matrix = slice_matrix.transpose()

    # calculate n: the number of q points, q_1 = qmin, q_n = qmax
    unique_q_vec = np.unique(iqmod.mod_q)
    minq_index = np.argmin(np.abs(qi_matrix[:, 0] - q_min))
    maxq_index = np.argmin(np.abs(qi_matrix[:, 0] - q_max))
    num_q = maxq_index - minq_index + 1
    print(f'Number of Q between {q_min} and {q_max} is {num_q}')

    # Unique wavelength and output b-vector
    unique_wavelength_vec = np.unique(iqmod.wavelength)
    b_factor_vec = np.zeros_like(unique_wavelenght_vec)

    # Loop over all wavelength
    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # select I(Q, wl=lambda_i)
        qi_matrix = wqi_matrix[wqi_matrix[:, 0] == lambda_i][:, 1:]

        # calculate n: the number of q points, q_1 = qmin, q_n = qmax
        # since I(Q, wavelength) are on mesh grid as a requirement on the input,
        # thus N is same for all wavelength
        np.testing.assert_allclose(unique_q_vec, qi_matrix[:, 0])

        # calculate b[i_wl] = -1 / num_q \sum_{min_q}^{max_q} [ I(q, ref_l) - I(q, wl) ]
        b_fact = 0.
        for i_q in range(minq_index, maxq_index + 1):
            b_fact += ref_wl_matrix[i_q, 2] - qi_matrix[i_q, 1]
        # b_factor = -1/num_q * b_fact 
        b_factor_vec[i_wl] = b_fact / (-num_q)
    # END-FOR

    return b_factor_vec


def select_min_incoherence(b_factor_vec, wavelength_vec):

    # find wavelength (lambda_i) has smallest b
    min_b_wl = blabla

    # choose minimum b wavelength as reference wavelength
    reference_wavelength_index = min_b_wl_index

    # recalculate reference wavelngth

    return


def correct_intensity_error(iqmod, b_factor_vec, ref_qie_matrix):

    # Retrieve values
    unique_wavelength_vec = np.unique(iqmod.wavelength)
    unique_q_vec = np.unique(iqmod.mod_q)
    num_unique_q = unique_q_vec.shape[0]
    num_points = iqmod.intensity.shape[0]
    assert unique_wavelength_vec.shape[0] * unique_q_vec.shape[0] == num_points

    # Create ND matrix for slicing
    # Create output
    if iqmod.delta_mod_q is None:
        complete_wqi_matrix = np.array([iqmod.wavelength, iqmod.mod_q, iqmod.intensity, iqmod.error])
        corrected_iqmod = IQmod(np.zeros_like(iqmod.intensity), np.zeros_like(iqmod.intensity),
                                np.zeros_like(iqmod.intensity), np.zeros_like(iqmod.intensity))
    else:
        complete_wqi_matrix = np.array([iqmod.wavelength, iqmod.mod_q, iqmod.intensity, iqmod.error, iqmod.delta_mod_q])
        corrected_iqmod = IQmod(np.zeros_like(iqmod.intensity), np.zeros_like(iqmod.intensity),
                                np.zeros_like(iqmod.intensity), np.zeros_like(iqmod.intensity), np.zeros_like(iqmod.intensity))
    complete_wqi_matrix = complete_wqi_matrix.transpose()

    # Loop over all the wavelength
    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # select I(q, wl=wl_i) and Error(q, wl)
        wqi_matrix = complete_wqi_matrix[complete_wqi_matrix[:, 0] == lambda_i]
        assert wqi_matrix.shape[0] == unique_q_vec.shape[0]

        # correct intensity I'(q, wl_i) = I(q, wl_i) - b(wl_i)
        corrected_intensity_vec = wqi_matrix[:, 2] - b_factor_vec[i_wl]

        # apply
        start_index = i_wl * num_unique_q
        stop_index = start_index + num_unique_q

        corrected_iqmod.wavelength[start_index:stop_index] = wqi_matrix[:, 0]
        corrected_iqmod.mod_q[start_index:stop_index] = wqi_matrix[:, 1]
        corrected_iqmod.intensity[start_index:stop_index] = corrected_intensity_vec
        if iqmod.delta_mod_q:
            corrected_iqmod.delta_mod_q[start_index:stop_index] = wqi_matrix[:, 4]

    return corrected_iqmod

