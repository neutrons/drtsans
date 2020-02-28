import numpy as np
import pytest
from drtsans.iq import determine_1d_linear_bins, determine_1d_log_bins
from tests.unit.new.drtsans.i_of_q_binning_tests_data import get_gold_2d_linear_bins, get_gold_1d_log_bins


# This module supports testing data for issue #239.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/263

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# Some tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data
def test_log_bins_backward_compatible():
    """Test log bins determination with 'old' API

    Method determine_1d_log_bins() has been refactored from its previous version by adding more
    method parameters.  While by default value, this method shall be backward compatible such that
    with x-min, x-max and step-per-decade defined, it shall generate a set of bins same as before.

    Here by using data '1D_bin_log_wedget_no_sub_no_wt' (from
    https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/5423db9b77dfd4911bf799c247530865/
    eqsans_tof_q_binning_tests_R5.xlsx), this method is tested

    Returns
    -------

    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    step_per_decade = 10  # 10 steps per decade

    log_bins = determine_1d_log_bins(q_min, q_max, step_per_decade, even_decade=True)
    gold_edges, gold_centers = get_gold_1d_log_bins()
    np.testing.assert_allclose(log_bins.edges, gold_edges, rtol=5.E-4)
    np.testing.assert_allclose(log_bins.centers, gold_centers, rtol=5.E-4)

    return


# Test EXCEL can be found at
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/5423db9b77dfd4911bf799c247530865/
# eqsans_tof_q_binning_tests_R5.xlsx
def test_linear_bin_determination():
    """Test linear bin determination from '2D_bin_no_sub_no_wt'

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


# Tests are from https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/263
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/ca6d7fe686595dda52ddef1f93d005f1/
# log_bin_definition_testsR1.xlsx
# All tests' gold data below have 3 columns as bin's left boundary, center and right boundary,
# which are exactly same in the test Excel file

# Example 1: Qmin and Qmax are even decade; Qmin and Qmax are bin centers
#            User defines number of bins per decade
gold_log_bins_example1 = np.array([
    [0.000089, 0.000100, 0.000113],
    [0.000113, 0.000126, 0.000142],
    [0.000142, 0.000158, 0.000179],
    [0.000179, 0.000200, 0.000225],
    [0.000225, 0.000251, 0.000284],
    [0.000284, 0.000316, 0.000357],
    [0.000357, 0.000398, 0.000450],
    [0.000450, 0.000501, 0.000566],
    [0.000566, 0.000631, 0.000713],
    [0.000713, 0.000794, 0.000897],
    [0.000897, 0.001000, 0.001129],
    [0.001129, 0.001259, 0.001422],
    [0.001422, 0.001585, 0.001790],
    [0.001790, 0.001995, 0.002254],
    [0.002254, 0.002512, 0.002837],
    [0.002837, 0.003162, 0.003572],
    [0.003572, 0.003981, 0.004496],
    [0.004496, 0.005012, 0.005661],
    [0.005661, 0.006310, 0.007126],
    [0.007126, 0.007943, 0.008972],
    [0.008972, 0.010000, 0.011295],
    [0.011295, 0.012589, 0.014219],
    [0.014219, 0.015849, 0.017901],
    [0.017901, 0.019953, 0.022536],
    [0.022536, 0.025119, 0.028371],
    [0.028371, 0.031623, 0.035717],
    [0.035717, 0.039811, 0.044965],
    [0.044965, 0.050119, 0.056607],
    [0.056607, 0.063096, 0.071264],
    [0.071264, 0.079433, 0.089716],
    [0.089716, 0.100000, 0.112202]])


# Example 2: Qmin and Qmax are even decade; Qmin and Qmax are on bin boundaries.
#            User defines number of bins per decade
gold_log_bins_example2 = np.array([
    [0.000100, 0.000112, 0.000127],
    [0.000127, 0.000141, 0.000160],
    [0.000160, 0.000178, 0.000201],
    [0.000201, 0.000224, 0.000253],
    [0.000253, 0.000282, 0.000318],
    [0.000318, 0.000355, 0.000401],
    [0.000401, 0.000447, 0.000505],
    [0.000505, 0.000562, 0.000635],
    [0.000635, 0.000708, 0.000800],
    [0.000800, 0.000891, 0.001007],
    [0.001007, 0.001122, 0.001267],
    [0.001267, 0.001413, 0.001595],
    [0.001595, 0.001778, 0.002009],
    [0.002009, 0.002239, 0.002529],
    [0.002529, 0.002818, 0.003183],
    [0.003183, 0.003548, 0.004007],
    [0.004007, 0.004467, 0.005045],
    [0.005045, 0.005623, 0.006351],
    [0.006351, 0.007079, 0.007996],
    [0.007996, 0.008913, 0.010066],
    [0.010066, 0.011220, 0.012673],
    [0.012673, 0.014125, 0.015954],
    [0.015954, 0.017783, 0.020085],
    [0.020085, 0.022387, 0.025286],
    [0.025286, 0.028184, 0.031833],
    [0.031833, 0.035481, 0.040075],
    [0.040075, 0.044668, 0.050451],
    [0.050451, 0.056234, 0.063514],
    [0.063514, 0.070795, 0.079960],
    [0.079960, 0.089125, 0.100000]])


# Example 3: Qmin and Qmax are exact as user specified (may not on decade).
#            User defines total number of bins.
gold_log_bins_example3 = np.array([
    [0.001500, 0.001582, 0.001671],
    [0.001671, 0.001759, 0.001858],
    [0.001858, 0.001957, 0.002066],
    [0.002066, 0.002176, 0.002298],
    [0.002298, 0.002420, 0.002556],
    [0.002556, 0.002692, 0.002843],
    [0.002843, 0.002993, 0.003161],
    [0.003161, 0.003329, 0.003516],
    [0.003516, 0.003703, 0.003910],
    [0.003910, 0.004118, 0.004349],
    [0.004349, 0.004580, 0.004837],
    [0.004837, 0.005093, 0.005379],
    [0.005379, 0.005665, 0.005982],
    [0.005982, 0.006300, 0.006653],
    [0.006653, 0.007007, 0.007399],
    [0.007399, 0.007792, 0.008229],
    [0.008229, 0.008666, 0.009152],
    [0.009152, 0.009638, 0.010179],
    [0.010179, 0.010719, 0.011320],
    [0.011320, 0.011922, 0.012590],
    [0.012590, 0.013259, 0.014002],
    [0.014002, 0.014746, 0.015573],
    [0.015573, 0.016400, 0.017319],
    [0.017319, 0.018239, 0.019262],
    [0.019262, 0.020285, 0.021422],
    [0.021422, 0.022560, 0.023825],
    [0.023825, 0.025090, 0.026497],
    [0.026497, 0.027904, 0.029469],
    [0.029469, 0.031033, 0.032774],
    [0.032774, 0.034514, 0.036398]])


# Example 4: Qmin and Qmax are exact as user specified (may not on decade).
#            User defines number of bins per decade
gold_log_bins_example4 = np.array([
    [0.000437, 0.000470, 0.000507],
    [0.000507, 0.000545, 0.000588],
    [0.000588, 0.000631, 0.000682],
    [0.000682, 0.000732, 0.000790],
    [0.000790, 0.000848, 0.000915],
    [0.000915, 0.000983, 0.001061],
    [0.001061, 0.001139, 0.001229],
    [0.001229, 0.001320, 0.001424],
    [0.001424, 0.001529, 0.001651],
    [0.001651, 0.001772, 0.001913],
    [0.001913, 0.002054, 0.002217],
    [0.002217, 0.002380, 0.002569],
    [0.002569, 0.002758, 0.002977],
    [0.002977, 0.003196, 0.003450],
    [0.003450, 0.003704, 0.003998],
    [0.003998, 0.004292, 0.004633],
    [0.004633, 0.004974, 0.005369],
    [0.005369, 0.005764, 0.006222],
    [0.006222, 0.006680, 0.007210],
    [0.007210, 0.007741, 0.008355],
    [0.008355, 0.008970, 0.009683],
    [0.009683, 0.010395, 0.011221],
    [0.011221, 0.012047, 0.013003],
    [0.013003, 0.013960, 0.015069],
    [0.015069, 0.016178, 0.017463],
    [0.017463, 0.018748, 0.020237],
    [0.020237, 0.021726, 0.023452],
    [0.023452, 0.025177, 0.027177],
    [0.027177, 0.029177, 0.031494],
    [0.031494, 0.033812, 0.036398]])


def test_example1():
    """Example 1 from ...

    10^c_min and 10^_max will be on decade.
    And 10^c_min and 10^c_max will be on bin boundaries but not bin centers (example 1)


    Returns
    -------

    """
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_per_decade = 10

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example1, q_max_example1, n_bins_per_decade, decade_on_center=True,
                                     even_decade=True)
    # verify bin center
    np.testing.assert_allclose(test_bin.centers, gold_log_bins_example1[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_bin.edges[:-1], gold_log_bins_example1[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_bin.edges[1:], gold_log_bins_example1[:, 2], rtol=1e-7, atol=1e-6)


def test_example2():
    """Example 2

    Example 2 has the same initial condition as example 1's.
    But the 10^c_min and 10^c_max will be on bin boundaries but not bin centers (example 1)

    Returns
    -------

    """
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_per_decade = 10

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example1, q_max_example1, n_bins_per_decade, decade_on_center=False,
                                     even_decade=True)
    # verify bin center
    np.testing.assert_allclose(test_bin.centers, gold_log_bins_example2[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_bin.edges[:-1], gold_log_bins_example2[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_bin.edges[1:], gold_log_bins_example2[:, 2], rtol=1e-7, atol=1e-6)


def test_example3():
    """Example 3

    X min and X max are exact as user specified (not on decade).  User defines total number of bins (30).
    Both Xmin and X max are on the bin boundaries.

    Test from Example 3 in log_bin_definition_testsR1.xlsx

    Returns
    -------

    """
    # Example 3
    q_min_example3 = 0.0015
    q_max_example3 = 0.036398139
    n_bins_example3 = 30

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bins = determine_1d_log_bins(q_min_example3, q_max_example3, n_bins=n_bins_example3,
                                      decade_on_center=False, even_decade=False)

    # verify bin center
    np.testing.assert_allclose(test_bins.centers, gold_log_bins_example3[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min (left boundary)
    np.testing.assert_allclose(test_bins.edges[:-1], gold_log_bins_example3[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max (right boundary)
    np.testing.assert_allclose(test_bins.edges[1:], gold_log_bins_example3[:, 2], rtol=1e-7, atol=1e-6)

    return


def test_example4():
    """Example 4

    X min and X max shall be on the bin boundaries;
    User specified X min and X max are used as X min/max exactly.

    Test from Example 4 in log_bin_definition_testsR1.xlsx

    Returns
    -------

    """
    # Example 4
    q_min_example4 = 0.0004367265
    q_max_example4 = 0.0363981393
    n_bins_per_decade = 10

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example4, q_max_example4, n_bins_per_decade, decade_on_center=False,
                                     even_decade=False)

    # verify bin center: center column of the gold data
    np.testing.assert_allclose(test_bin.centers, gold_log_bins_example4[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min: left column of the gold data
    np.testing.assert_allclose(test_bin.edges[:-1], gold_log_bins_example4[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max: right column of the gold data
    np.testing.assert_allclose(test_bin.edges[1:], gold_log_bins_example4[:, 2], rtol=1e-7, atol=1e-6)

    return


if __name__ == '__main__':
    pytest.main([__file__])
