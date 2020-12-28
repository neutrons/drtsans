import pytest
from drtsans.dataobjects import IQmod

from drtsans.tof.eqsans.elastic_reference_normalization import determine_common_mod_q_range
from drtsans.tof.eqsans.elastic_reference_normalization import normalize_by_elastic_reference, determine_reference_wavelength_q1d
from drtsans.tof.eqsans.elastic_reference_normalization import normalize_by_elastic_reference, determine_reference_wavelength_q1d

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
    vec_intensity, vec_error, vec_q, vec_wl = create_testing_iq1d(num_q, num_wl)
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
    vec_intensity, vec_error, vec_q, vec_wl = create_testing_iq1d(num_q, num_wl)
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
    assert test_ref_lambda_vec.shape == (num_q, 2)

    # Verify the values
    assert pytest.approx(test_ref_lambda_vec[0][0],  100.)
    assert test_ref_lambda_vec[0][1] == 0.4
    assert pytest.approx(test_ref_lambda_vec[1][1], 0.3)
    assert test_ref_lambda_vec[2][1] == 0.2
    for i in range(3, num_q):
        assert test_ref_lambda_vec[i][1] == 0.1, f'{i}-th Q has a wrong reference wavelength ' \
                                                 f'{test_ref_lambda_vec[i]}'


def create_testing_iq1d(num_q, num_wl):
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


if __name__ == '__main__':
    pytest.main([__file__])
