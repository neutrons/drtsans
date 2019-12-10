import numpy as np
from drtsans.dataobjects import IQazimuthal
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/iq.py
from drtsans.iq import determine_1d_linear_bins, _do_1d_no_weight_binning,\
    BinningMethod, BinningParams, bin_annular_into_q1d
from tests.unit.new.drtsans.i_of_q_binning_tests_data import generate_test_data, get_gold_theta_bins,\
    get_gold_azimuthal_values
import pytest


# This module supports testing data for issue #246.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/246

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data


def test_1d_annular_no_wt_prototype():
    """Test '1D_annular_no_sub_no_wt' by verifying prototype of 1D annular binning method

    Returns
    -------

    """
    theta_min = 0
    theta_max = 360.
    num_bins = 10

    q_min = 0.003
    q_max = 0.006

    # Generate testing data: Get Q2D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Calculate theta array
    theta_array = np.arctan2(qy_array, qx_array) * 180. / np.pi
    # convert -0 to -180 to 180 to 360
    theta_array[np.where(theta_array < 0)] += 360.

    # Calculate Q from Qx and Qy
    q_array = np.sqrt(qx_array**2 + qy_array**2)
    # calculate dQ from dQx and dQy
    dq_array = np.sqrt(dqx_array**2 + dqy_array**2)
    # Filter by q_min and q_max
    allowed_q_index = (q_array > q_min) & (q_array < q_max)

    # Generate theta bins
    theta_bins = determine_1d_linear_bins(theta_min, theta_max, num_bins)

    # binning
    binned_iq = _do_1d_no_weight_binning(theta_array[allowed_q_index],
                                         dq_array[allowed_q_index],
                                         intensities[allowed_q_index],
                                         sigmas[allowed_q_index],
                                         theta_bins.centers, theta_bins.edges)

    # Check bins
    gold_theta_edges, gold_theta_centers = get_gold_theta_bins()
    np.testing.assert_allclose(theta_bins.centers, gold_theta_centers, rtol=1.e-5)
    np.testing.assert_allclose(theta_bins.edges, gold_theta_edges, rtol=1.e-5)

    # Check theta (azimuthal angle)
    # assert abs(theta_array[0] - 141.026949) < 5E-2, 'Azimuthal angle check'
    gold_theta_array = get_gold_azimuthal_values()
    num_test_data = gold_theta_array.shape[0]

    np.testing.assert_allclose(theta_array[:num_test_data], gold_theta_array, rtol=6.e-4, err_msg='Azimuthal vectors')

    # Verify I(Q), sigma I(Q) and dQ
    assert binned_iq.intensity[1] == pytest.approx(63.66666667, abs=1E-8), 'Binned intensity is wrong'
    assert binned_iq.error[1] == pytest.approx(3.257470048, abs=1E-8), 'Binned sigma I is wrong'
    assert binned_iq.delta_mod_q[1] == pytest.approx(1.154E-02, abs=1E-5), 'Binned Q resolution {} ' \
                                                                           'is incorrect comparing to {}.' \
                                                                           ''.format(binned_iq.delta_mod_q[1],
                                                                                     0.01154)


def test_1d_annular_no_wt():
    """Test '1D_annular_no_sub_no_wt'

    Returns
    -------

    """
    # Initialize range of theta angle and Q
    theta_min = 0
    theta_max = 360.
    num_bins = 10

    q_min = 0.003
    q_max = 0.006

    # Generate testing data: Get Q2D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test the high level method
    # Define input data
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)

    # Annular binning
    theta_binning = BinningParams(theta_min, theta_max, num_bins)
    binned_iq = bin_annular_into_q1d(test_i_q, theta_binning, q_min, q_max, BinningMethod.NOWEIGHT)

    # Verify I(Q), sigma I(Q) and dQ
    assert binned_iq.intensity[1] == pytest.approx(63.66666667, abs=1E-8), 'Binned intensity is wrong'
    assert binned_iq.error[1] == pytest.approx(3.257470048, abs=1E-8), 'Binned sigma I is wrong'
    assert binned_iq.delta_mod_q[1] == pytest.approx(1.154E-02, abs=1E-5), 'Binned Q resolution {} ' \
                                                                           'is incorrect comparing to {}.' \
                                                                           ''.format(binned_iq.delta_mod_q[1],
                                                                                     0.01154)


if __name__ == '__main__':
    pytest.main([__file__])
