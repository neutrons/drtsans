# Test drtsans.tof.eqsans.incoherence_correction_1d
import pytest
from drtsans.dataobjects import IQmod
from drtsans.tof.eqsans.elastic_reference_normalization import reshape_q_wavelength_matrix
from drtsans.tof.eqsans.elastic_reference_normalization import determine_common_mod_q_range_mesh
from drtsans.tof.eqsans.elastic_reference_normalization import determine_reference_wavelength_q1d_mesh
import numpy as np


def test_incoherence_inelastic_correction():
    """
    (Serving as a prototype)
    """
    select_min_incoh = False

    # generate testing data
    test_iq1d = generate_test_data()

    # convert to mesh grid I(Q) and delta I(Q)
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_iq1d)

    # determine q min and q ma that exists in all I(q, wl) ...
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    b_array, ref_wl_ie = calculate_b_factors(wl_vec, q_vec, i_array, error_array, select_min_incoh,
                                             qmin_index, qmax_index)
    # # determine the reference wavelength
    # ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
    #                                                     qmin_index, qmax_index)
    #
    # # calculate b(lambda_i)
    # b_array = calculate_b_vector_prototype(wl_vec, q_vec, i_array, qmin_index, qmax_index, ref_wl_ie,
    #                                        calculate_delta_b=not select_min_incoh)
    # # If JSON parameter “selectMinIncoh” is true
    # if select_min_incoh:
    #     # it seems that there is no test case for this
    #     ref_wl_ie = reselect_reference_wl(b_vector)
    #     b_array = calculate_b_vector_prototype(wl_vec, q_vec, i_array, qmin_index, qmax_index, ref_wl_ie,
    #                                            calculate_delta_b=True)
    #     # verify
    #     assert b_array[np.isnan(b_array)] >= 0.
    # else:
    #     # set minimum wavelength bin to nan
    #     if abs(b_array[0, 0]) < 1E-10:
    #         b_array[0, 0] = np.nan
    #         b_array[1, 0] = np.nan

    # check
    np.testing.assert_allclose(b_array[0], generate_expected_b_factors(), verbose=True)

    # Do correction
    corrected = correct_intensity_error_prototype(i_array, error_array, b_array[0])

    # verify
    corrected_i_array = corrected[0].flatten()
    np.testing.assert_allclose(corrected_i_array, generate_expected_corrected_intensities(),
                               verbose=True)


def test_incoherence_inelastic_correction_select_min_incoherence():
    select_min_incoh = True

    # generate testing data
    test_iq1d = generate_data_select_min_incoherence()

    # convert to mesh grid I(Q) and delta I(Q)
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_iq1d)

    # determine q min and q ma that exists in all I(q, wl) ...
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    b_array, ref_wl_ie = calculate_b_factors(wl_vec, q_vec, i_array, error_array, select_min_incoh,
                                             qmin_index, qmax_index)

    print(b_array)
    print(ref_wl_ie)


def correct_intensity_error_prototype(wavelength_vec, q_vec, intensity_array, error_array, b_vector,
                                      qmin_index, qmax_index, ref_wl_ie):

    # Correct intensity: I'(q, wl) = I(q, wl) - b(wl)
    for i_wl in range(b_vector.shape[0]):
        intensity_array[:, i_wl] -= b_vector[i_wl]

    # Correct error
    # correction on error cannot be in place
    corrected_error2_array = np.zeros_like(error_array)

    for i_wl in range(len(wavelength_vec)):

        n_q = qmax_index + 1 - qmin_index

        for q_j in range(len(q_vec)):

            # delta I'(j, i) = (delta I(j, i))^2 * (1 - 2/N) + ...
            error2_j = error_array[q_j, i_wl]**2 * (1 - 2./n_q)

            for q_k in range(qmin_index, qmax_index+1):
                # delta I'(j, i) += 1/N^2 [(Ref_delta_I(k))^2 + (delta I(k, i)^2]: i for wl
                error2_j += 1./n_q**2 * (ref_wl_ie.error_vec[q_k]**2 + error_array[q_k, i_wl]**2)

            corrected_error2_array[q_j, i_wl] = error2_j

    return intensity_array, np.sqrt(corrected_error2_array)


def calculate_b_vector_prototype(wl_vec, q_vec, i_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                 calculate_delta_b):

    # row 0: b(wl)  row1: delta b(wl)
    b_array = np.ndarray(shape=(2, len(wl_vec)), dtype='float')

    for i_wl in range(len(wl_vec)):

        n_q = qmax_index - qmin_index + 1

        b_i = 0.
        delta_b_i = 0.

        for j_q in range(qmin_index, qmax_index+1):
            # b = -1/N sum_{q_j} (refI(q_j) - I(q_j, wl)
            b_i += ref_wl_ie.intensity_vec[j_q] - i_array[j_q, i_wl]
            # delta b^2 = 1/N^2 sum_{q_j} [(ref_delta_I(q_j)^2 + (delta I(q_j, wl)^2]
            delta_b_i += 1. / n_q**2 * (ref_wl_ie.error_vec[j_q]**2 + error_array[j_q, i_wl]**2)

        b_array[0, i_wl] = - b_i / n_q




    # # set minimum wavelength bin to nan
    # if abs(b_vector[0]) < 1E-10:
    #     b_vector[0] = np.nan

    return b_array


def generate_test_data():
    # Generate test data given in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/
    # /b3b4038f44443385afe4252bb2316be3/inelastic_incoherent_avg_example.xlsx
    # Denoted as TEST1

    # Intensity vector
    intensity_vec = np.array([
        0.1, 0.13, np.nan, np.nan, np.nan,  # q = 0.01
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])

    # Q vector
    vec_q = np.arange(1, 21) * 0.01
    vec_q = np.repeat(vec_q, 5)
    print(f'q: {vec_q.shape}')

    # Wavelength vector
    wavelength_vec = np.arange(1, 6) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)
    print(f'lambda: {wavelength_vec.shape}')

    # Error
    error_vec = np.sqrt(intensity_vec)

    # Construct IQmod
    i_of_q = IQmod(intensity=intensity_vec,
                   error=error_vec,
                   mod_q=vec_q,
                   wavelength=wavelength_vec)

    return i_of_q


def generate_expected_b_factors():
    # Generate test data given in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/
    # /b3b4038f44443385afe4252bb2316be3/inelastic_incoherent_avg_example.xlsx
    # Denoted as TEST1

    # Expected B vectors
    b_factor_vec = np.array([np.nan, -0.03, -0.05, -0.04, -0.01])

    return b_factor_vec


def generate_expected_corrected_intensities():

    # Expected corrected intensities
    corrected_intensity_vec = np.array([
        0.1, 0.1, np.nan, np.nan, np.nan,  # q = 0.01
        0.1, 0.1, 0.1, np.nan, np.nan,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        np.nan, 0.1, 0.1, 0.1, 0.1,
        np.nan, 0.1, 0.1, 0.1, 0.1,
        np.nan, np.nan, 0.1, 0.1, 0.1,
        np.nan, np.nan, 0.1, 0.1, 0.1,
        np.nan, np.nan, np.nan, 0.1, 0.1,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])

    return corrected_intensity_vec


def generate_data_select_min_incoherence():
    """Generate data to test the methods to calculate B(ref wavelength)
    with 'select min incoherence' option

    Returns
    -------

    """
    intensity_vec = np.array([
        0.1, 0.13, np.nan, np.nan, np.nan,  # q = 0.01
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.09, 0.11,
        0.1, 0.13, 0.15, 0.04, 0.11,
        0.1, 0.13, 0.15, 0.04, 0.11,
        0.1, 0.13, 0.15, 0.04, 0.11,
        0.1, 0.13, 0.15, 0.04, 0.11,
        np.nan, 0.13, 0.15, 0.04, 0.11,
        np.nan, 0.13, 0.15, 0.04, 0.11,
        np.nan, np.nan, 0.15, 0.04, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])

    # Q vector
    vec_q = np.arange(1, 21) * 0.01
    vec_q = np.repeat(vec_q, 5)
    print(f'q: {vec_q.shape}')

    # Wavelength vector
    wavelength_vec = np.arange(1, 6) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)
    print(f'lambda: {wavelength_vec.shape}')

    # Error
    error_vec = np.sqrt(intensity_vec)

    # Construct IQmod
    i_of_q = IQmod(intensity=intensity_vec,
                   error=error_vec,
                   mod_q=vec_q,
                   wavelength=wavelength_vec)

    return i_of_q


def calculate_b_factors(wl_vec, q_vec, i_array, error_array, select_min_incoherence, qmin_index, qmax_index):

    # determine reference wavelength and calculate b vectors

    # determine the reference wavelength
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                        qmin_index, qmax_index, 0)

    # calculate b(lambda_i) and delta b(lambda_i) if it is final
    b_array = calculate_b_vector_prototype(wl_vec, q_vec, i_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                           calculate_delta_b=not select_min_incoherence)

    # If JSON parameter “selectMinIncoh” is true
    if select_min_incoherence and np.argmin(b_array[0]) > 0:
        # reselect reference wavelength to the wavelength bin with minimum b
        ref_wl_index = np.argmin(b_array[0])
        print(f'[DEBUG...] re-selected reference wavelength index = {ref_wl_index}')
        # (re)determine the reference wavelengths' intensities and errors
        ref_wl_ie = determine_reference_wavelength_q1d_mesh(ref_wl_index, wl_vec, q_vec, i_array, error_array,
                                                            qmin_index, qmax_index)
        # (re)calculate b array
        b_array = calculate_b_vector_prototype(wl_vec, q_vec, i_array, qmin_index, qmax_index, ref_wl_ie,
                                               calculate_delta_b=True)
        # verify
        assert b_array[np.isnan(b_array)] >= 0.

    return b_array, ref_wl_ie


if __name__ == '__main__':
    pytest.main([__file__])

