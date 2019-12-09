import numpy as np
import pytest


# Equation 11.31
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


def define_log_bin(x_min=0.0001, x_max=0.1, n_bins=10, decade_on_center=False):
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
    c_min = np.log10(x_min)  # c_min may not be on decade!
    c_max = np.ceil(np.log10(x_max))

    print('[TEST] C min = {}, C max = {}'.format(c_min, c_max))

    # Calculate total number of bins
    total_num_bins = (c_max - c_min) * n_bins
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
    else:
        # x_min and 10^{c_max} on the first and last bin boundary
        # then first and last bin edges/boundaries are defined as
        # Q_min and Q_max (or X_min and X_max in generalized form)
        bin_edges[0] = x_min
        bin_edges[-1] = x_max
    # END-IF-ELSE

    return c_min, c_max, total_num_bins, delta_l, bin_centers, bin_edges


def test_example1():
    """

    Returns
    -------

    """
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_example1 = 10

    # Example 1: Max/Min on Bin Centers'
    print('[TEST] Example 1: Max/Min on Bin Centers')
    test_set = define_log_bin(q_min_example1, q_max_example1, n_bins_example1, True)
    print(test_set[4])
    print(gold_log_bins_example1[:, 1])
    diffs = np.abs(test_set[4] - gold_log_bins_example1[:, 1])
    print('Max Mean Diff = {}'.format(np.max(diffs)))

    print(test_set[5])
    print(gold_log_bins_example1[:, 0])

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


def next_test_example2():
    # Test data for both Example 1 and Example 2
    q_min_example1 = 0.0001
    q_max_example1 = 0.036398139
    n_bins_example1 = 10
    # Example 2
    print('[TEST] Example 2: Max/Min on Bin Boundaries')
    test_set = define_log_bin(q_min_example1, q_max_example1, n_bins_example1, False)
    # Verify with expected value
    gold_c_min = -4
    gold_c_max = -1
    gold_n_bins = 30
    gold_delta_l = 0.1
    assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-10, '{} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, '{} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-10, '{} != {}'.format(test_set[3], gold_delta_l)
    # verify bin center
    np.testing.assert_allclose(test_set[4], gold_log_bins_example2[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example2[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    set_a = test_set[5][1:]
    set_b = gold_log_bins_example2[:, 2]
    print(np.abs(set_a - set_b))
    print(np.max(np.abs(set_a - set_b)))
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example2[:, 2], rtol=1e-7, atol=1e-6)


def next_test_example3():
    # Example 3
    print('[TEST] Example 3: Min Q and Max Q on bin boundaries')
    q_min_example3 = 0.0015
    q_max_example3 = 0.036398139
    n_bins_example3 = 10

    test_set = define_log_bin(q_min_example3, q_max_example3, n_bins_example3, False)
    # Verify with expected value
    gold_c_min = -2.823908741
    gold_c_max = -1.438920817
    gold_n_bins = 30
    gold_delta_l = 0.046166264
    assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-10, '{} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, '{} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-10, '{} != {}'.format(test_set[3], gold_delta_l)
    # verify bin center
    np.testing.assert_allclose(test_set[4], gold_log_bins_example3[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example3[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example3[:, 2], rtol=1e-7, atol=1e-6)

    # Example 4: Cmin and Cmax are not on decades!
    q_min_example4 = 0.000436726
    q_max_example4 = 0.036398139

    test_set = define_log_bin(q_min_example4, q_max_example4, n_bins_example3, False)
    # Verify with expected value
    gold_c_min = -3.359790455
    gold_c_max = -1.438920817
    gold_n_bins = 30
    gold_delta_l = 0.064028988
    assert abs(test_set[0] - gold_c_min) < 1E-10, '{} != {}'.format(test_set[0], gold_c_min)
    assert abs(test_set[1] - gold_c_max) < 1E-10, '{} != {}'.format(test_set[1], gold_c_max)
    assert abs(test_set[2] - gold_n_bins) < 1E-10, '{} != {}'.format(test_set[2], gold_n_bins)
    assert abs(test_set[3] - gold_delta_l) < 1E-10, '{} != {}'.format(test_set[3], gold_delta_l)
    # verify bin center
    np.testing.assert_allclose(test_set[4], gold_log_bins_example4[:, 1], rtol=1e-7, atol=1e-6)
    # verify bin boundaries min
    np.testing.assert_allclose(test_set[5][:-1], gold_log_bins_example4[:, 0], rtol=1e-7, atol=1e-6)
    # verify bin boundaries max
    np.testing.assert_allclose(test_set[5][1:], gold_log_bins_example4[:, 2], rtol=1e-7, atol=1e-6)

    return


if __name__ == '__main__':
    pytest.main([__file__])