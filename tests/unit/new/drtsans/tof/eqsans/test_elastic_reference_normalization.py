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
    # Get data:
    test_i_of_q = create_testing_iq1d()

    # Calculate Qmin and Qmax
    q_min, q_max = determine_common_mod_q_range(test_i_of_q)
    assert q_min == 0.03
    assert q_max == 0.13

    # Calculate reference wavelength
    ref_q_wl_vec = determine_reference_wavelength_q1d(test_i_of_q)
    # vector Q
    vec_q = ref_q_wl_vec[:, 0]
    np.testing.assert_allclose(vec_q, np.arange(1, 21) * 0.01)
    # wavelength
    vec_wl = ref_q_wl_vec[:, 1]
    np.testing.assert_allclose(vec_wl[:13], 1.)
    np.testing.assert_allclose(vec_wl[13:], np.array([2., 2., 3., 3., 4., 4., 5.]))

    # Calculate scale factor
    k_vec, delta_k_vec, p_vec, s_vec, unique_wl_vec, ref_q_wl_vec = calculate_scale_factor(test_i_of_q, q_min, q_max)

    # Verify scale factor
    assert k_vec.shape == (5,)
    assert delta_k_vec.shape == (5,)

    # Normalize
    # TODO / FIXME - shall not be here
    normalize_intensity(test_i_of_q, k_vec, delta_k_vec, ref_q_wl_vec, p_vec, s_vec,
                        unique_wl_vec)


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
        0.1110, 0.13, 0.15, 0.14, 0.11,
        0.1111, 0.13211, 0.15, 0.14, 0.11,
        np.nan, 0.13212, 0.15, 0.14, 0.11,
        np.nan, 0.13213, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])
    print(f'intensity: {intensity_vec.shape}')

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


if __name__ == '__main__':
    pytest.main([__file__])
