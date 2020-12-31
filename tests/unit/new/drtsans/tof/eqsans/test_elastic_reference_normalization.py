import pytest
from drtsans.dataobjects import IQmod

from drtsans.tof.eqsans.elastic_reference_normalization import determine_common_mod_q_range
from drtsans.tof.eqsans.elastic_reference_normalization import determine_reference_wavelength_q1d
from drtsans.tof.eqsans.elastic_reference_normalization import calculate_scale_factor
from drtsans.tof.eqsans.elastic_reference_normalization import normalize_intensity

import numpy as np


def test_determine_q_range():
    """Test the method to determine the common q-range among all the wavelength vectors

    Test data
        - L wavelength, 1., 2., 3.,
        - N Q on regular grids for each wavelength
        - For each wavelength, there are P (P <= N) consecutive arbitrary finite values
        WL/Q   q1  q2  q3  ... ... q_N
        wl1    nan v12 v13 ... ... nan
        wl2    nan nan v23 ... ....
        ...
        wl_L
    """
    # Create test data set and manually set the range of valid wavelength for each Q
    num_wl = 10
    num_q = 100
    vec_intensity, vec_error, vec_q, vec_wl = create_configurable_testing_iq1d(num_q, num_wl)
    # set wl index = 1 with min q = 300, max q = 70000
    vec_intensity[num_q:num_q + 2] = np.nan
    vec_intensity[num_q * 2 - 3:num_q * 2] = np.nan
    # create the IQmod
    test_i_q1d = IQmod(intensity=vec_intensity, error=vec_error,
                       mod_q=vec_q, wavelength=vec_wl)

    min_q, max_q = determine_common_mod_q_range(test_i_q1d)

    assert min_q == 300.
    assert max_q == 9700.


def test_determine_reference_wavelength_q1d():
    """Test method to determine the reference wavelength for each Q in I(Q, wavelength)

    Test data
        - L wavelength, 1., 2., 3.,
        - N Q on regular grids for each wavelength
        - For each wavelength, there are P (P <= N) consecutive arbitrary finite values
        WL/Q   q1  q2  q3  ... ... q_N
        wl1    nan nan nan ... ... nan
        wl2    nan nan v13 ... ....
        wl3    nan v32 v33
        wl4    v41 v42 v43 ... ..
        wl_L
    """
    # Create test data set and manually set the range of valid wavelength for each Q
    num_wl = 10
    num_q = 100
    vec_intensity, vec_error, vec_q, vec_wl = create_configurable_testing_iq1d(num_q, num_wl)
    # set minimum wavelength range: first 3 from wl_i, such as (wl_1, wl_2, wl_3)
    # q1, q2, q3
    vec_intensity[0:3] = np.nan
    vec_intensity[num_q:num_q+2] = np.nan
    vec_intensity[num_q*2:num_q*2+1] = np.nan

    # create the IQmod
    test_i_q1d = IQmod(intensity=vec_intensity, error=vec_error,
                       mod_q=vec_q, wavelength=vec_wl)

    # Call method to calculate reference
    test_ref_lambda_vec = determine_reference_wavelength_q1d(test_i_q1d)

    # Test: shape... (num Q, 2) as Q and wavelength
    assert test_ref_lambda_vec.shape == (num_q, 4)

    # Verify the values
    assert pytest.approx(test_ref_lambda_vec[0][0],  100.)
    assert test_ref_lambda_vec[0][1] == 0.4
    assert pytest.approx(test_ref_lambda_vec[1][1], 0.3)
    assert test_ref_lambda_vec[2][1] == 0.2
    for i in range(3, num_q):
        assert test_ref_lambda_vec[i][1] == 0.1, f'{i}-th Q has a wrong reference wavelength ' \
                                                 f'{test_ref_lambda_vec[i]}'


def test_calculate_scale_factor():
    """Test method to calculate scale factor K(lambda) and error delta K(lambda)

    Test data: using data from attached Excel
    """
    # Get testing data and gold data
    test_i_of_q, gold_k_vec, gold_intensity_vec, gold_error_vec = create_testing_iq1d()

    # Calculate Qmin and Qmax
    q_min, q_max = determine_common_mod_q_range(test_i_of_q)
    assert q_min == 0.07
    assert q_max == 0.11

    # Calculate reference wavelength
    ref_q_wl_vec = determine_reference_wavelength_q1d(test_i_of_q)
    # vector Q
    vec_q = ref_q_wl_vec[:, 0]
    np.testing.assert_allclose(vec_q, np.arange(1, 21) * 0.01)
    # wavelength
    vec_wl = ref_q_wl_vec[:, 1]
    np.testing.assert_allclose(vec_wl[6:], 3.)
    np.testing.assert_allclose(vec_wl[0:6], np.array([6., 6., 5., 5., 4., 4.]))

    # Calculate scale factor
    k_vec, delta_k_vec, p_vec, s_vec, unique_wl_vec, ref_q_wl_vec = calculate_scale_factor(test_i_of_q, q_min, q_max)

    # Verify scale factor, P[3] and S[3]
    # considering numerical precision between original Excel with more effective decimals than the copied version
    # to python text
    assert k_vec.shape == (4,)
    assert delta_k_vec.shape == (4,)
    np.testing.assert_allclose(k_vec, gold_k_vec, atol=1.E-3)
    assert p_vec[3] == pytest.approx(1266.145, 0.005)
    assert s_vec[3] == pytest.approx(1139.531, 0.005)
    assert delta_k_vec[3] == pytest.approx(0.2063066, abs=0.00002)  # numerical precision difference

    # Normalize
    # Due to the numerical precision issue when import data from Excel
    # it is more precise to use expected k vector
    corrected_iqmod = normalize_intensity(test_i_of_q, gold_k_vec,
                                          ref_q_wl_vec, p_vec, s_vec,
                                          unique_wl_vec, q_min, q_max)

    # verify: shape
    assert corrected_iqmod
    assert corrected_iqmod.mod_q.shape == test_i_of_q.mod_q.shape

    # verify: mod_q and wavelength
    # original: each row has same wavelength; output: each column has same wavelength
    # solution: transpose original
    np.testing.assert_allclose(corrected_iqmod.wavelength,
                               test_i_of_q.wavelength.reshape((20, 4)).transpose().flatten())
    # transpose as above
    np.testing.assert_allclose(corrected_iqmod.mod_q,
                               test_i_of_q.mod_q.reshape((20, 4)).transpose().flatten())

    # verify: corrected intensity
    # gold = gold_intensity_vec.reshape((20, 4)).transpose().flatten()
    # for i in range(80):
    #     buf = f'{corrected_iqmod.intensity[i]} -  {gold[i]}'
    #     if corrected_iqmod.intensity[i] is not None and gold[i] is not None:
    #         buf += f'  = {corrected_iqmod.intensity[i] - gold[i]}'
    #     print(buf)

    np.testing.assert_allclose(corrected_iqmod.intensity,
                               gold_intensity_vec.reshape((20, 4)).transpose().flatten(),
                               equal_nan=True,
                               atol=1E-3)

    # verify: corrected intensity error
    gold_errors = gold_error_vec.reshape((20, 4)).transpose().flatten()
    for i in range(80):
        buf = f'{corrected_iqmod.error[i]} -  {gold_errors[i]}'
        if corrected_iqmod.intensity[i] is not None and gold_errors[i] is not None:
            buf += f'  = {corrected_iqmod.error[i] - gold_errors[i]}'
        print(buf)

    np.testing.assert_allclose(corrected_iqmod.error,
                               gold_errors,
                               equal_nan=True,
                               atol=1E-5)


def create_configurable_testing_iq1d(num_q, num_wl):
    # create a testing I(Q)
    vec_q = (np.arange(num_q).astype('float') + 1.) * 100.
    vec_wl = (np.arange(num_wl).astype('float') + 1.) * 0.1

    # make a full list of Q and lambda such in order as
    # q1, wl1
    # q2, wl1
    # ...
    # q1, wl2
    vec_q = np.tile(vec_q, num_wl)
    vec_wl = np.repeat(vec_wl, num_q)

    # intensity
    vec_intensity = np.zeros(shape=(num_q * num_wl,), dtype='float') + 123.
    vec_error = np.sqrt(vec_intensity)

    return vec_intensity, vec_error, vec_q, vec_wl


def create_testing_iq1d():
    """Create a test data I(Q, wavelength) as the attached EXCEL spreadsheet attached in gitlab story
    Returns
    -------

    """
    # Intensity vector
    intensity_vec = np.array([
        np.nan, np.nan, np.nan, 743.802,
        np.nan, np.nan, np.nan, 459.184,
        np.nan, np.nan, 332.410, 249.307,
        np.nan, np.nan, 177.515, 133.136,
        np.nan, 89.796, 97.959, 73.469,
        np.nan, 51.985, 56.711, 42.533,
        28.727, 31.600, 34.473, 25.855,
        18.262, 20.088, 21.914, 16.435,
        12.076, 13.283, 14.491, 10.868,
        8.264, 9.091, 9.917, 7.438,
        5.827, 6.410, 6.993, 5.244,
        4.217, 4.638, 5.060, np.nan,
        3.121, 3.433, 3.745, np.nan,
        2.356, 2.592, 2.828, np.nan,
        1.811, 1.992, 2.173, np.nan,
        1.413, 1.555, np.nan, np.nan,
        1.119, 1.230, np.nan, np.nan,
        0.896, np.nan, np.nan, np.nan,
        0.727, np.nan, np.nan, np.nan,
        0.595, np.nan, np.nan, np.nan,
    ])
    print(f'intensity: {intensity_vec.shape}')

    # Error
    error_vec = np.array([
        np.nan, np.nan, np.nan, 27.273,
        np.nan, np.nan, np.nan, 21.429,
        np.nan, np.nan, 18.232, 15.789,
        np.nan, np.nan, 13.323, 11.538,
        np.nan, 16.792, 9.897, 8.571,
        np.nan, 10.611, 7.531, 6.522,
        5.360, 4.605, 5.871, 5.085,
        4.273, 4.374, 4.681, 4.054,
        3.475, 3.640, 3.807, 3.297,
        2.875, 2.985, 3.149, 2.727,
        2.414, 2.470, 2.644, 2.290,
        2.053, 2.095, 2.249, np.nan,
        1.767, 1.772, 1.935, np.nan,
        1.535, 1.522, 1.682, np.nan,
        1.346, 1.322, 1.474, np.nan,
        1.189, 1.161, np.nan, np.nan,
        1.058, 1.028, np.nan, np.nan,
        0.947, np.nan, np.nan, np.nan,
        0.852, np.nan, np.nan, np.nan,
        0.771, np.nan, np.nan, np.nan,
    ])

    # Q vector
    vec_q = np.arange(1, 21) * 0.01
    vec_q = np.repeat(vec_q, 4)
    print(f'q: {vec_q.shape}')

    # Wavelength vector
    wavelength_vec = np.arange(3, 7) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)
    print(f'lambda: {wavelength_vec.shape}')

    # Construct IQmod
    i_of_q = IQmod(intensity=intensity_vec,
                   error=error_vec,
                   mod_q=vec_q,
                   wavelength=wavelength_vec)

    # Expected K vector
    expected_k_vec = np.array([1.0000, 0.9090909090909, 0.833333333333333, 1.11111111111111])

    expected_corrected_intensities = np.array([
        np.nan, np.nan, np.nan, 826.446,
        np.nan, np.nan, np.nan, 510.204,
        np.nan, np.nan, 277.008, 277.008,
        np.nan, np.nan, 147.929, 147.929,
        np.nan, 81.633, 81.633, 81.633,
        np.nan, 47.259, 47.259, 47.259,
        28.727, 28.727, 28.727, 28.727,
        18.262, 18.262, 18.262, 18.262,
        12.076, 12.076, 12.076, 12.076,
        8.264, 8.264, 8.264, 8.264,
        5.827, 5.827, 5.827, 5.827,
        4.217, 4.217, 4.217, np.nan,
        3.121, 3.121, 3.121, np.nan,
        2.356, 2.356, 2.356, np.nan,
        1.811, 1.811, 1.811, np.nan,
        1.413, 1.413, np.nan, np.nan,
        1.119, 1.119, np.nan, np.nan,
        0.896, np.nan, np.nan, np.nan,
        0.727, np.nan, np.nan, np.nan,
        0.595, np.nan, np.nan, np.nan,
    ])

    expected_corrected_errors = np.array([
        np.nan, np.nan, np.nan, 156.415,
        np.nan, np.nan, np.nan, 97.679,
        np.nan, np.nan, 50.281, 54.344,
        np.nan, np.nan, 27.900, 30.312,
        np.nan, 9.476, 16.357, 17.901,
        np.nan, 7.210, 10.308, 11.380,
        5.360, 5.621, 4.534, 4.788,
        4.273, 4.482, 4.241, 4.708,
        3.475, 3.645, 3.513, 3.958,
        2.875, 3.015, 2.875, 3.263,
        2.414, 2.532, 2.374, 2.708,
        2.053, 2.154, 2.011, np.nan,
        1.767, 1.853, 1.701, np.nan,
        1.535, 1.610, 1.459, np.nan,
        1.346, 1.411, 1.268, np.nan,
        1.189, 1.247, np.nan, np.nan,
        1.058, 1.109, np.nan, np.nan,
        0.947, np.nan, np.nan, np.nan,
        0.852, np.nan, np.nan, np.nan,
        0.771, np.nan, np.nan, np.nan,
    ])

    return i_of_q, expected_k_vec, expected_corrected_intensities, expected_corrected_errors


if __name__ == '__main__':
    pytest.main([__file__])
