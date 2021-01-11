# Test drtsans.tof.eqsans.incoherence_correction_1d
import pytest
from drtsans.dataobjects import IQmod
from drtsans.tof.eqsans.elastic_reference_normalization import reshape_q_wavelength_matrix
from drtsans.tof.eqsans.elastic_reference_normalization import determine_common_mod_q_range_mesh
from drtsans.tof.eqsans.elastic_reference_normalization import determine_reference_wavelength_q1d_mesh
from drtsans.tof.eqsans.incoherence_correction_1d import calculate_b_error_b, calculate_b_factors
from drtsans.tof.eqsans.incoherence_correction_1d import correct_intensity_error, correct_incoherence_inelastic_1d
import numpy as np


def test_calculate_b_factor():
    """Test calculate B factor and reference wavelengths by Changwoo's test case
    """
    select_min_incoh = False

    # generate testing data
    test_iq1d = generate_test_data()

    # convert to mesh grid I(Q) and delta I(Q)
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_iq1d)

    # determine q min and q ma that exists in all I(q, wl) ...
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # test unit method
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                        qmin_index, qmax_index, 0)
    test_b_array = calculate_b_error_b(wl_vec, i_array, error_array, qmin_index, qmax_index,
                                       ref_wl_ie, True)
    # verify
    np.testing.assert_allclose(test_b_array[0], generate_expected_b_factors())

    # call prototype
    b_array, ref_wl_ie = calculate_b_factors_prototype(wl_vec, q_vec, i_array, error_array, select_min_incoh,
                                                       qmin_index, qmax_index)

    # verify
    np.testing.assert_allclose(b_array[0], generate_expected_b_factors())
    # test delta B by 2 various implementation
    np.testing.assert_allclose(test_b_array[1], b_array[1])


def test_calculate_b_factor_select_min_incoherence():
    """Test calculate B factor and reference wavelength with option 'select minimum incoherence'
    """
    select_min_incoh = True

    # generate testing data
    test_iq1d = generate_data_select_min_incoherence()

    # convert to mesh grid I(Q) and delta I(Q)
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_iq1d)

    # determine q min and q ma that exists in all I(q, wl) ...
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # calculate B factor and error with select minimum incoherence flag on
    b_array, ref_wl_ie = calculate_b_factors(wl_vec, q_vec, i_array, error_array, select_min_incoh,
                                             qmin_index, qmax_index)

    # verify
    assert np.argmin(b_array[0]) == 3
    assert b_array[0, 3] == 0.

    # call prototype to compare with
    b_array_prototype, ref_wl_ie_p = calculate_b_factors_prototype(wl_vec, q_vec, i_array, error_array,
                                                                   select_min_incoh, qmin_index, qmax_index)
    np.testing.assert_allclose(b_array, b_array_prototype)


def test_incoherence_inelastic_correction():
    """Test methods to correct I(Q1D) accouting wavelength dependent incoherence
    inelastic scattering
    """
    select_min_incoh = False

    # generate testing data
    test_iq1d = generate_test_data()

    # convert to mesh grid I(Q) and delta I(Q)
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_iq1d)

    # determine q min and q ma that exists in all I(q, wl) ...
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # calculate B factors and errors
    b_array, ref_wl_ie = calculate_b_factors(wl_vec, q_vec, i_array, error_array, select_min_incoh,
                                             qmin_index, qmax_index)
    # verify
    np.testing.assert_allclose(b_array[0], generate_expected_b_factors(), verbose=True)

    # correct intensities and errors
    corrected_intensities, corrected_errors = correct_intensity_error(wl_vec, q_vec, i_array, error_array,
                                                                      b_array, qmin_index, qmax_index,
                                                                      ref_wl_ie)
    # Verify intensities
    corrected_i_array = corrected_intensities.flatten()

    np.testing.assert_allclose(corrected_i_array, generate_expected_corrected_intensities(),
                               verbose=True)

    # Do correction by prototypes
    corrected = correct_intensity_error_prototype(wl_vec, q_vec, i_array, error_array, b_array[0],
                                                  qmin_index, qmax_index, ref_wl_ie)

    # Compare results from correct_intensity_error and prototypes
    np.testing.assert_allclose(corrected_errors, corrected[1])

    # Test overall workflow
    corrected_i_of_q = correct_incoherence_inelastic_1d(test_iq1d, False)
    np.testing.assert_allclose(corrected_i_of_q.intensity, generate_expected_corrected_intensities())


def correct_intensity_error_prototype(wavelength_vec, q_vec, intensity_array, error_array, b_vector,
                                      qmin_index, qmax_index, ref_wl_ie):

    # Correct intensity: I'(q, wl) = I(q, wl) - b(wl)
    for i_wl in range(b_vector.shape[0]):
        intensity_array[:, i_wl] -= b_vector[i_wl]

    # Correct error
    # correction on error cannot be in place
    corrected_error2_array = np.zeros_like(error_array)
    n_q = qmax_index + 1 - qmin_index

    for i_wl in range(len(wavelength_vec)):
        for q_j in range(len(q_vec)):

            if qmin_index <= q_j <= qmax_index:
                # q_j is within q1...qN:
                # delta I'(j, i) = (delta I(j, i))^2 * (1 - 2/N) + ...
                error2_j = error_array[q_j, i_wl]**2 * (1 - 2./n_q)

                for q_k in range(qmin_index, qmax_index+1):
                    # delta I'(j, i) += 1/N^2 [(Ref_delta_I(k))^2 + (delta I(k, i)^2]: i for wl
                    error2_j += 1./n_q**2 * (ref_wl_ie.error_vec[q_k]**2 + error_array[q_k, i_wl]**2)
            else:
                # q_j is outside q1...qN
                # [delta I'(q_j, wl)]^2 = [delta I(q_j, wl)]^2 + ...
                error2_j = error_array[q_j, i_wl]**2
                # [delta I'(q_j, wl)]^2 += 1/N^2 sum_{q_k in qmin:qmax+1}[(RefDelta I(q_k))^2 + (delta I(q_k, wl)^2]
                for q_k in range(qmin_index, qmax_index+1):
                    error2_j += 1. / n_q**2 * (ref_wl_ie.error_vec[q_k]**2 + error_array[q_k, i_wl]**2)

            corrected_error2_array[q_j, i_wl] = error2_j

    return intensity_array, np.sqrt(corrected_error2_array)


def calculate_b_vector_prototype(wl_vec, i_array, error_array, qmin_index, qmax_index, ref_wl_ie,
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

            if calculate_delta_b:
                # delta b^2 = 1/N^2 sum_{q_j} [(ref_delta_I(q_j)^2 + (delta I(q_j, wl)^2]
                delta_b_i += 1. / n_q**2 * (ref_wl_ie.error_vec[j_q]**2 + error_array[j_q, i_wl]**2)

        b_array[0, i_wl] = - b_i / n_q
        if calculate_delta_b:
            b_array[1, i_wl] = np.sqrt(delta_b_i)

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

    # Wavelength vector
    wavelength_vec = np.arange(1, 6) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)

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
    b_factor_vec = np.array([0., 0.03, 0.05, 0.04, 0.01])

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
        np.nan, np.nan, np.nan, 0.1, 0.1,
        np.nan, np.nan, np.nan, np.nan, 0.1,
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

    # Wavelength vector
    wavelength_vec = np.arange(1, 6) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)

    # Error
    error_vec = np.sqrt(intensity_vec)

    # Construct IQmod
    i_of_q = IQmod(intensity=intensity_vec,
                   error=error_vec,
                   mod_q=vec_q,
                   wavelength=wavelength_vec)

    return i_of_q


def calculate_b_factors_prototype(wl_vec, q_vec, i_array, error_array, select_min_incoherence, qmin_index, qmax_index):

    # determine reference wavelength and calculate b vectors

    # determine the reference wavelength
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                        qmin_index, qmax_index, 0)

    # calculate b(lambda_i) and delta b(lambda_i) if it is final
    b_array = calculate_b_vector_prototype(wl_vec, i_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                           calculate_delta_b=not select_min_incoherence)

    # If JSON parameter “selectMinIncoh” is true
    if select_min_incoherence and np.argmin(b_array[0]) > 0:
        # reselect reference wavelength to the wavelength bin with minimum b
        ref_wl_index = np.argmin(b_array[0])
        # (re)determine the reference wavelengths' intensities and errors
        ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                            qmin_index, qmax_index, ref_wl_index)
        # (re)calculate b array
        b_array = calculate_b_vector_prototype(wl_vec, i_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                               calculate_delta_b=True)
        # verify
        assert b_array[np.isfinite(b_array)].min() >= -1E-20, f'B array has negative values: {b_array}'

    return b_array, ref_wl_ie


if __name__ == '__main__':
    pytest.main([__file__])
