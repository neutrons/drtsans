import numpy as np
import pytest
from drtsans.iq import determine_1d_linear_bins, determine_1d_log_bins
from tests.unit.new.drtsans.i_of_q_binning_tests_data import get_gold_2d_linear_bins, get_gold_1d_log_bins


# This module supports testing data for issue #239.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/263

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# Some tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data
def test_log_bins_from_wedge_no_wt():
    """Test generating log bins from '1D_bin_log_wedget_no_sub_no_wt'
    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    step_per_decade = 10  # 10 steps per decade

    log_bins = determine_1d_log_bins(q_min, q_max, step_per_decade)
    gold_edges, gold_centers = get_gold_1d_log_bins()
    np.testing.assert_allclose(log_bins.edges, gold_edges, rtol=5.E-4)
    np.testing.assert_allclose(log_bins.centers, gold_centers, rtol=5.E-4)

    return


# Test EXCEL can be found at
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/5423db9b77dfd4911bf799c247530865/
# eqsans_tof_q_binning_tests_R4.xlsx
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
# Example 1: Equation 11.31
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


# Example 2
# Equation 11.29 - 11.30
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


# Example 3
# Equation 11.29 - 11.30
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


# Example 4
# Equation 11.29 - 11.30
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


def determine_log_bin_prototype(x_min=0.0001, x_max=0.1, n_bins=10, decade_on_center=False, even_decade=True):
    """Determine logarithm bins

    Including bin edge and bin center

    Parameters
    ----------
    x_min : float
        Positive float for minimum
    x_max : float
        Positive float
    n_bins : int
        Positive integer for number of step per decade
    decade_on_center : bool
        Flag to have the min X and max X on bin center; Otherwise, they will be on bin boundary

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        bin centers, bin edges
    """
    # Calculate C min and max on decade and contain X min and X max in the range
    if decade_on_center or even_decade:
        c_min = np.floor(np.log10(x_min))  # c_min may not be on decade!
        c_max = np.ceil(np.log10(x_max))
    else:
        c_min = np.log10(x_min)  # c_min may not be on decade!
        c_max = np.log10(x_max)

    print('[TEST] C min = {}, C max = {}'.format(c_min, c_max))

    # Calculate total number of bins
    total_num_bins = (int(np.ceil(c_max) - np.floor(c_min))) * n_bins
    print('[TEST] Number of total bins = {}'.format(total_num_bins))

    # Calculate Delta L
    delta_l = (c_max - c_min) / total_num_bins
    print('[TEST] Bin size DeltaL = {}'.format(delta_l))
    # Define an array of k, i.e., [0, 1, 2, ...]
    if decade_on_center:
        # x_min and 10^{c_max} on the first and last bin center
        # Equation 11.31: number of bins will be N + 1 for bin on the center
        total_num_bins += 1
        vec_k = np.arange(total_num_bins).astype(float)
        bin_centers = np.power(10, delta_l * vec_k + c_min)
    else:
        # x_min and 10^{c_max} on the first and last bin boundary
        # Equation 11.29
        vec_k = np.arange(total_num_bins).astype(float)
        bin_centers = np.power(10, delta_l * (vec_k + 0.5) + c_min)
    # END-IF-ELSE

    # Calculate bin boundaries from bin center
    # Equation 11.30
    bin_edges = np.ndarray(shape=(bin_centers.shape[0] + 1,), dtype=float)
    bin_edges[1:-1] = 0.5 * (bin_centers[:-1] + bin_centers[1:])

    if decade_on_center:
        # x_min and 10^{c_max} are on the first and last bin center
        # then first and last bin edges/boundaries are defined as
        # 10^{C_min - deltaL / 2 } and 10^{C_max + deltaL / 2}
        # according to the paragraph after equation 11.31
        bin_edges[0] = np.power(10, c_min - 0.5 * delta_l)
        bin_edges[-1] = np.power(10, c_max + 0.5 * delta_l)
    elif even_decade:
        # x_min and 10^{c_max} on the first and last bin boundary
        # then first and last bin edges/boundaries are defined as
        # Q_min and Q_max (or X_min and X_max in generalized form)
        bin_edges[0] = 10**c_min
        bin_edges[-1] = 10**c_max
    else:
        # use user x min and x max
        bin_edges[0] = x_min
        bin_edges[-1] = x_max
    # END-IF-ELSE

    return c_min, c_max, total_num_bins, delta_l, bin_centers, bin_edges


def test_example1():
    """Example 1 from ...

    10^c_min and 10^_max will be on decade.
    And 10^c_min and 10^c_max will be on bin boundaries but nto bin centers (example 1)


    Returns
    -------

    """
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_example1 = 10

    # Test with prototype: will be removed after all test cases are passed
    test_set = determine_log_bin_prototype(q_min_example1, q_max_example1, n_bins_example1, True)
    # Verify with expected value
    gold_c_min = -4
    gold_c_max = -1
    gold_n_bins = 31
    gold_delta_l = 0.1
    assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-10, '{} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, '{} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-10, '{} != {}'.format(test_set[3], gold_delta_l)
    # verify bin center
    np.testing.assert_allclose(test_set[4], gold_log_bins_example1[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example1[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example1[:, 2], rtol=1e-7, atol=1e-6)

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example1, q_max_example1, n_bins_example1, decade_on_center=True,
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
    But the 10^c_min and 10^c_max will be on bin boundaries but nto bin centers (example 1)

    Returns
    -------

    """
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_example1 = 10

    # TODO - prototype will be deleted after all tess are passed
    # prototype test
    test_set = determine_log_bin_prototype(q_min_example1, q_max_example1, n_bins_example1, decade_on_center=False,
                                           even_decade=True)
    # Verify with expected value
    gold_c_min = -4
    gold_c_max = -1
    gold_n_bins = 30
    gold_delta_l = 0.1
    assert abs(test_set[0] - gold_c_min) < 1E-10, 'cmin {} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-10, 'cmax {} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, 'number of bins {} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-10, 'delta L {} != {}'.format(test_set[3], gold_delta_l)

    # verify bin center: same number ....
    np.testing.assert_allclose(test_set[4], gold_log_bins_example2[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example2[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example2[:, 2], rtol=1e-7, atol=1e-6)

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example1, q_max_example1, n_bins_example1, decade_on_center=False,
                                     even_decade=True)
    # verify bin center
    np.testing.assert_allclose(test_bin.centers, gold_log_bins_example2[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_bin.edges[:-1], gold_log_bins_example2[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_bin.edges[1:], gold_log_bins_example2[:, 2], rtol=1e-7, atol=1e-6)


def test_example3():
    """Test for example 3 from log_bin_definition_testsR1.xlsx

    In this case, total number of bins are not calculated from bin per decade but given by user

    Returns
    -------

    """
    # Example 3
    print('[TEST] Example 3: Min Q and Max Q on bin boundaries')
    q_min_example3 = 0.0015
    q_max_example3 = 0.036398139
    n_bins_example3 = 30

    test_bins = determine_1d_log_bins(q_min_example3, q_max_example3, n_bins=None, decade_on_center=False,
                                      even_decade=False)

    # verify bin center
    np.testing.assert_allclose(test_bins.centers, gold_log_bins_example4[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min (left boundary)
    np.testing.assert_allclose(test_bins.edges[:-1], gold_log_bins_example4[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max (right boundary)
    np.testing.assert_allclose(test_bins.edges[1:], gold_log_bins_example4[:, 2], rtol=1e-7, atol=1e-6)

    # test_set = determine_log_bin_prototype(q_min_example3, q_max_example3,  , decade_on_center=False,
    #                                        even_decade=False)
    # # Verify with expected value
    # gold_c_min = -2.823908741
    # gold_c_max = -1.438920817
    # gold_n_bins = 30
    # gold_delta_l = 0.046166264
    # assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    # assert abs(test_set[1] - gold_c_max) < 1E-7, 'CMAX {} != {}'.format(test_set[1], gold_c_max)
    # assert abs(test_set[2] - gold_n_bins) < 1E-10, 'TOTAL BINS {} != {}'.format(test_set[2], gold_n_bins)
    # assert abs(test_set[3] - gold_delta_l) < 1E-10, 'DELTA L {} != {}'.format(test_set[3], gold_delta_l)
    # # verify bin center
    # np.testing.assert_allclose(test_set[4], gold_log_bins_example3[:, 1], rtol=1e-7, atol=1e-6)
    # # verify bin boundaries min
    # np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example3[:, 0], rtol=1e-7, atol=1e-6)
    # # verify bin boundaries max
    # np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example3[:, 2], rtol=1e-7, atol=1e-6)
    #
    # # Example 4: Cmin and Cmax are not on decades!
    # q_min_example4 = 0.000436726
    # q_max_example4 = 0.036398139
    #
    # test_set = determine_log_bin_prototype(q_min_example4, q_max_example4, n_bins_example3, False)
    # # Verify with expected value
    # gold_c_min = -3.359790455
    # gold_c_max = -1.438920817
    # gold_n_bins = 30
    # gold_delta_l = 0.064028988
    # assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    # assert abs(test_set[1] - gold_c_max) < 1E-10, '{} != {}'.format(test_set[1], gold_c_max)
    # assert abs(test_set[2] - gold_n_bins) < 1E-10, '{} != {}'.format(test_set[2], gold_n_bins)
    # assert abs(test_set[3] - gold_delta_l) < 1E-10, '{} != {}'.format(test_set[3], gold_delta_l)
    # # verify bin center
    # np.testing.assert_allclose(test_set[4], gold_log_bins_example4[:, 1], rtol=1e-7, atol=1e-6)
    # # verify bin boundaries min
    # np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example4[:, 0], rtol=1e-7, atol=1e-6)
    # # verify bin boundaries max
    # np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example4[:, 2], rtol=1e-7, atol=1e-6)

    return


def test_example4():
    """Test case from example 4

    X min and X max shall be on the bin boundaries;
    User specified X min and X max are used as X min/max exactly.

    Returns
    -------

    """
    # Example 4
    q_min_example4 = 0.0004367265
    q_max_example4 = 0.0363981393
    n_bins_example = 10

    # Prototype test: TODO - delete after all tests are passed
    test_set = determine_log_bin_prototype(q_min_example4, q_max_example4, n_bins_example, decade_on_center=False,
                                           even_decade=False)
    # Verify with expected value
    gold_c_min = -3.359790455
    gold_c_max = -1.438920817
    gold_n_bins = 30
    gold_delta_l = 0.064028988
    assert abs(test_set[0] - gold_c_min) < 1E-7, 'CMIN {} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-7, 'CMAX {} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, 'TOTAL BINS {} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-7, 'DELTA L {} != {}'.format(test_set[3], gold_delta_l)
    # verify bin center
    np.testing.assert_allclose(test_set[4], gold_log_bins_example4[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example4[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example4[:, 2], rtol=1e-7, atol=1e-6)

    # Test drtsans.determine_bins.determine_1d_log_bins
    test_bin = determine_1d_log_bins(q_min_example4, q_max_example4, n_bins_example, decade_on_center=False,
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
