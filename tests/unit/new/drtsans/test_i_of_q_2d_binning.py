import numpy as np
from drtsans.dataobjects import IQazimuthal
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/iq.py
from drtsans.iq import determine_1d_linear_bins, BinningMethod, _do_2d_weighted_binning,\
    _do_2d_no_weight_binning, bin_intensity_into_q2d
from tests.unit.new.drtsans.i_of_q_binning_tests_data import generate_test_data, get_gold_2d_linear_bins
import pytest

# This module supports testing data for issue #239.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/245

# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data


def test_2d_linear_bin():
    """Test '2D_bin_no_sub_no_wt' and '2D_bin_no_sub_wt'

    2D linear bin no sub pixel with weighted and no-weight summation

    Returns
    -------
    None

    """
    # Calculate and determine the bin edges
    # range of binned (Qx, Qy) is taken from William's Excel
    qx_min = -0.007573828
    qx_max = 0.006825091
    qy_min = -0.005051412
    qy_max = 0.00607504

    qx_bins = determine_1d_linear_bins(qx_min, qx_max, 5)
    qy_bins = determine_1d_linear_bins(qy_min, qy_max, 5)

    # verify
    gold_x_centers, gold_y_centers = get_gold_2d_linear_bins()

    np.testing.assert_allclose(qx_bins.centers, gold_x_centers, atol=5E-6)
    np.testing.assert_allclose(qy_bins.centers, gold_y_centers, atol=5E-6)

    # Check X
    assert qx_bins.edges[1] == pytest.approx(-0.004694044, abs=1E-8)
    assert qx_bins.edges[2] == pytest.approx(-0.001814261, abs=1E-8)
    # Check Y
    assert qy_bins.edges[1] == pytest.approx(-0.002826, abs=1E-6)
    assert qy_bins.edges[2] == pytest.approx(-0.000601, abs=1E-6)

    # Bin 2D No-weight
    # Get Q1D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test for no-weight binning
    binned_iq_2d = _do_2d_no_weight_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                            qx_bins.edges, qy_bins.edges)

    # Verify Qx and Qy
    assert qx_bins.centers[1] == pytest.approx(-0.003254, abs=1.E-6), 'Qx is not correct'
    assert qy_bins.centers[1] == pytest.approx(-0.001713, abs=1.E-6), 'Qy is not correct'

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    assert binned_iq_2d[0][1][1] == pytest.approx(67., abs=1E-6), 'I(Qx, Qy) is incorrect'
    assert binned_iq_2d[1][1][1] == pytest.approx(4.725815626, abs=1E-8), 'sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    # correct: 3.2999999999999996e-05
    assert binned_iq_2d[2][1][1] == pytest.approx(0.00816, abs=1E-5), 'dQx {} is incorrect comparing to {}.' \
                                                                      ''.format(binned_iq_2d[2][1][1], 0.00816)
    assert binned_iq_2d[3][1][1] == pytest.approx(0.00816, abs=1E-5), 'dQy {}is incorrect comparing to {}.' \
                                                                      ''.format(binned_iq_2d[3][1][1], 0.00816)

    # Bin 2D Weighted
    # Test for weighted-binning
    binned_iq_2d = _do_2d_weighted_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                           qx_bins.edges, qy_bins.edges)

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    # test value: 56.86602493293357
    assert binned_iq_2d[0][1][1] == pytest.approx(56.8660, abs=1E-4), 'Weighted-binned I(Qx, Qy) is incorrect'
    assert binned_iq_2d[1][1][1] == pytest.approx(4.353773265, abs=1E-8), \
        'Weighted-binned sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    assert binned_iq_2d[2][1][1] == pytest.approx(0.00816, abs=1E-5), 'dQx is incorrect'
    # correct: 1.71877860186208e-05
    assert binned_iq_2d[3][1][1] == pytest.approx(0.00816, abs=1E-5), 'dQy is incorrect'

    # Test API for high level method
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)
    binned_iq_2d_wt = bin_intensity_into_q2d(test_i_q, qx_bins, qy_bins, BinningMethod.WEIGHTED)
    # verify
    np.testing.assert_allclose(binned_iq_2d[0], binned_iq_2d_wt.intensity, atol=1E-10)

    return


if __name__ == '__main__':
    pytest.main([__file__])
