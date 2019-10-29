import numpy as np
from drtsans.iq import determine_1d_linear_bins, determine_1d_log_bins, do_1d_no_weight_binning,\
    bin_iq_into_logarithm_q1d, BinningMethod, do_2d_weighted_binning, do_2d_no_weight_binning

# This test implements issue #169 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/tree/169_bin_q1d
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

intensities_matrix = np.array([[[93, 60, 89, 32, 97],
                                [43, 61, 82, 97, 55],
                                [78, 34, 50, 54, 67],
                                [98, 88, 37, 92, 97],
                                [72, 97, 100, 71, 39]],
                               [[76, 39, 51, 70, 61],
                                [64, 54, 78, 35, 30],
                                [67, 98, 100, 56, 79],
                                [97, 35, 41, 90, 45],
                                [30, 41, 68, 34, 51]],
                               [[78, 36, 46, 75, 91],
                                [64, 56, 92, 73, 60],
                                [74, 72, 69, 84, 87],
                                [36, 78, 40, 68, 72],
                                [59, 40, 39, 34, 85]]])

uncertainties_matrix = np.sqrt(intensities_matrix)

# William's input Qx
qx_matrix = np.array([[[-0.006134, -0.003254, -0.000374, 0.002505, 0.005385],
                       [-0.006134, -0.003254, -0.000374, 0.002505, 0.005385],
                       [-0.006134, -0.003254, -0.000374, 0.002505, 0.005385],
                       [-0.006134, -0.003254, -0.000374, 0.002505, 0.005385],
                       [-0.006134, -0.003254, -0.000374, 0.002505, 0.005385]],
                      [[-0.005936, -0.003149, -0.000362, 0.002425, 0.005211],
                       [-0.005936, -0.003149, -0.000362, 0.002425, 0.005211],
                       [-0.005936, -0.003149, -0.000362, 0.002425, 0.005211],
                       [-0.005936, -0.003149, -0.000362, 0.002425, 0.005211],
                       [-0.005936, -0.003149, -0.000362, 0.002425, 0.005211]],
                      [[-0.005751, -0.003051, -0.000351, 0.002349, 0.005049],
                       [-0.005751, -0.003051, -0.000351, 0.002349, 0.005049],
                       [-0.005751, -0.003051, -0.000351, 0.002349, 0.005049],
                       [-0.005751, -0.003051, -0.000351, 0.002349, 0.005049],
                       [-0.005751, -0.003051, -0.000351, 0.002349, 0.005049]]])

# Williams's input Qy
qy_matrix = np.array([[[0.004962, 0.004962, 0.004962, 0.004962, 0.004962],
                       [0.002737, 0.002737, 0.002737, 0.002737, 0.002737],
                       [0.000512, 0.000512, 0.000512, 0.000512, 0.000512],
                       [-0.001713, -0.001713, -0.001713, -0.001713, -0.001713],
                       [-0.003939, -0.003939, -0.003939, -0.003939, -0.003939]],
                      [[0.004802, 0.004802, 0.004802, 0.004802, 0.004802],
                       [0.002649, 0.002649, 0.002649, 0.002649, 0.002649],
                       [0.000495, 0.000495, 0.000495, 0.000495, 0.000495],
                       [-0.001658, -0.001658, -0.001658, -0.001658, -0.001658],
                       [-0.003812, -0.003812, -0.003812, -0.003812, -0.003812]],
                      [[0.004652, 0.004652, 0.004652, 0.004652, 0.004652],
                       [0.002566, 0.002566, 0.002566, 0.002566, 0.002566],
                       [0.000480, 0.000480, 0.000480, 0.000480, 0.000480],
                       [-0.001606, -0.001606, -0.001606, -0.001606, -0.001606],
                       [-0.003693, -0.003693, -0.003693, -0.003693, -0.003693]]])


dqx_matrix = np.array([[[0.000066, 0.000035, 0.000005, 0.000027, 0.000058],
                        [0.000066, 0.000035, 0.000005, 0.000027, 0.000058],
                        [0.000066, 0.000035, 0.000005, 0.000027, 0.000058],
                        [0.000066, 0.000035, 0.000005, 0.000027, 0.000058],
                        [0.000066, 0.000035, 0.000005, 0.000027, 0.000058]],
                       [[0.000062, 0.000033, 0.000004, 0.000025, 0.000055],
                        [0.000062, 0.000033, 0.000004, 0.000025, 0.000055],
                        [0.000062, 0.000033, 0.000004, 0.000025, 0.000055],
                        [0.000062, 0.000033, 0.000004, 0.000025, 0.000055],
                        [0.000062, 0.000033, 0.000004, 0.000025, 0.000055]],
                       [[0.000058, 0.000031, 0.000004, 0.000024, 0.000051],
                        [0.000058, 0.000031, 0.000004, 0.000024, 0.000051],
                        [0.000058, 0.000031, 0.000004, 0.000024, 0.000051],
                        [0.000058, 0.000031, 0.000004, 0.000024, 0.000051],
                        [0.000058, 0.000031, 0.000004, 0.000024, 0.000051]]])

dqy_matrix = np.array([[[0.000054, 0.000054, 0.000054, 0.000054, 0.000054],
                        [0.000030, 0.000030, 0.000030, 0.000030, 0.000030],
                        [0.000006, 0.000006, 0.000006, 0.000006, 0.000006],
                        [0.000019, 0.000019, 0.000019, 0.000019, 0.000019],
                        [0.000043, 0.000043, 0.000043, 0.000043, 0.000043]],
                       [[0.000050, 0.000050, 0.000050, 0.000050, 0.000050],
                        [0.000028, 0.000028, 0.000028, 0.000028, 0.000028],
                        [0.000006, 0.000006, 0.000006, 0.000006, 0.000006],
                        [0.000017, 0.000017, 0.000017, 0.000017, 0.000017],
                        [0.000040, 0.000040, 0.000040, 0.000040, 0.000040]],
                       [[0.000047, 0.000047, 0.000047, 0.000047, 0.000047],
                        [0.000026, 0.000026, 0.000026, 0.000026, 0.000026],
                        [0.000005, 0.000005, 0.000005, 0.000005, 0.000005],
                        [0.000016, 0.000016, 0.000016, 0.000016, 0.000016],
                        [0.000038, 0.000038, 0.000038, 0.000038, 0.000038]]])


def generate_test_data(q_dimension, drt_standard):
    """Generate test data

    Test data including I(Q), Q and dQ depending on 1D or 2D

    Parameters
    ----------
    q_dimension : int
        Scalar Q or Qx, Qy
    drt_standard: bool
        flag to convert test 3D data (detector view + wave length) to drt standard 1D

    Returns
    -------

    """
    # Check input
    if q_dimension not in [1, 2]:
        raise RuntimeError('Q-dimension must be 1 or 2')

    if q_dimension == 1:
        # Calculate scalar Q and dQ
        scalar_q_matrix = np.sqrt(qx_matrix ** 2 + qy_matrix ** 2)
        scalar_dq_matrix = np.sqrt(dqx_matrix ** 2 + dqy_matrix ** 2)
    else:
        # No-op
        scalar_dq_matrix = scalar_q_matrix = None

    # Define a None returning object
    returns = None
    if drt_standard:
        # convert intensities
        intensity_array = intensities_matrix.flatten()
        sigma_array = uncertainties_matrix.flatten()
        if q_dimension == 1:
            # 1D: scalar
            scalar_q_array = scalar_q_matrix.flatten()
            scalar_dq_array = scalar_dq_matrix.flatten()
            returns = intensity_array, sigma_array, scalar_q_array, scalar_dq_array
        else:
            qx_array = qx_matrix.flatten()
            qy_array = qy_matrix.flatten()
            dqx_array = dqx_matrix.flatten()
            dqy_array = dqy_matrix.flatten()
            returns = intensity_array, sigma_array, qx_array, dqx_array, qy_array, dqy_array
    else:
        # Raw matrix format
        if q_dimension == 1:
            # 1D: scalar Q
            returns = intensities_matrix, uncertainties_matrix, scalar_q_matrix, scalar_q_matrix
        elif q_dimension == 2:
            # 2D: Qx, Qy
            returns = intensities_matrix, uncertainties_matrix, qx_matrix, dqx_matrix, qy_matrix, dqy_matrix

    return returns


def get_gold_1d_linear_bins():
    """Get the gold array for 1D linear bins

    This is to test the method to create linear bins

    Returns
    -------
    ndarray, ndarray
        bin edges, bin centers
    """
    # Create bin edges (N + 1)
    edge_array = np.array([
        0.0000, 0.0010, 0.0020, 0.0030, 0.0040,
        0.0050, 0.0060, 0.0070, 0.0080, 0.0090, 0.0100])

    center_array = np.array([
        .0005, 0.0015, 0.0025, 0.0035, 0.0045, 0.0055, 0.0065,
        0.0075, 0.0085, 0.0095])

    return edge_array, center_array


def get_gold_1d_log_bins():
    """Get the gold array for 1D logarithm bins

    This is to test the method to create logarithm bins

    Returns
    -------
    ndarray, ndarray
        bin edges, bin centers
    """
    edge_array = np.array(
        [0.001000, 0.001267, 0.001595, 0.002009, 0.002529, 0.003183, 0.004007,
         0.005045, 0.006351, 0.007996, 0.010000])

    center_array = np.array(
        [0.001122, 0.001413, 0.001778, 0.002239, 0.002818, 0.003548, 0.004467,
         0.005623, 0.007079, 0.008913])

    return edge_array, center_array


def get_gold_2d_linear_bins():
    """Get the gold (Qx, Qy) bins

    Data are from William's Excel tests

    Returns
    -------
    ndafray, ndarray
        Qx centers, Qy centers
    """
    qx_center = np.array([-0.006134, -0.003254, -0.000374, 0.002505, 0.005385])

    qy_center = np.array([0.004962, 0.002737, 0.000512, -0.001713, -0.003939])

    return qx_center, qy_center


def test_1d_bin_linear_no_wt():
    """Test '1D_bin_linear_no_sub_wt'

    Test methods for 1D linear no-weight binning

    Returns
    -------
    None

    """
    q_min = 0.000
    q_max = 0.010
    num_bins = 10

    # Verify bin edges and bin center
    bin_centers, bin_edges = determine_1d_linear_bins(q_min, q_max, num_bins)
    gold_edges, gold_centers = get_gold_1d_linear_bins()

    assert np.allclose(bin_edges, gold_edges, 1.E-12)
    assert np.allclose(bin_centers, gold_centers, 1.E-12)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = do_1d_no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                        bin_centers, bin_edges)

    # Calculate and verify
    # I(0.0035) = 68.92857:    drtsans: 68.92857142857143
    assert abs(binned_iq.i[3] - 68.92857) < 2.E-6, 'I wrong'
    # di(0.0035)		2.218889
    assert abs(binned_iq.sigma[3] - 2.218889) < 2.E-6, 'sigma I wrong'
    # sigma_Q(0.0035) = 3.722E-05: This is off as it is the value from EXCEL with some error
    assert abs(binned_iq.dq[3] - 3.722E-05) < 2.E-5, 'Q resolution wrong'

    return


def test_1d_bin_log_no_wt():
    """Test '1D_bin_log_no_sub_no_wt'

    Test binning methods for 1D no-weight binning with log bins

    Returns
    -------

    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    num_steps_per_10 = 10  # 10 steps per decade

    # Verify bin edges and bin center
    bin_centers, bin_edges = determine_1d_log_bins(q_min, q_max, num_steps_per_10)
    gold_edges, gold_centers = get_gold_1d_log_bins()

    assert np.allclose(bin_edges, gold_edges, 5.E-4)
    assert np.allclose(bin_centers, gold_centers, 5.E-4)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = do_1d_no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                        bin_centers, bin_edges)

    # Verify: 2 I(Q) in bin: Q(3, 2, 3.1), Q(3, 2, 3.2)
    # I(0.0022) = 70.00000
    assert abs(binned_iq.i[3] - 70.00000) < 1.E-12, 'I wrong'
    # dI(0.0022) = 5.9160797831
    assert abs(binned_iq.sigma[3] - 5.9160797831) < 1.E-12, 'sigma I wrong'
    # sigma_Q(0.0022) = 2.529E-05: this value is from EXCEL with error in Q resolution
    # corrected value shall be   2.5112610804313703e-05
    assert abs(binned_iq.dq[3] - 2.529E-05) < 2.E-7, 'Q resolution wrong'

    # Test the high level method
    binned_iq = bin_iq_into_logarithm_q1d(intensities, sigmas, scalar_q_array, scalar_dq_array,
                                          num_steps_per_10, q_min, q_max, BinningMethod.NOWEIGHT)
    # I(0.0022) = 70.00000
    assert abs(binned_iq.i[3] - 70.00000) < 1.E-12, 'I wrong'

    return


def test_2d_linear_bin_no_wt():
    """Test '2D_bin_no_sub_no_wt'

    2D linear bin no sub pixel no weighing summation

    Returns
    -------

    """
    # Calculate and determine the bin edges
    # range of binned (Qx, Qy) is taken from William's Excel
    qx_min = -0.007573828
    qx_max = 0.006825091
    qy_min = -0.005051412
    qy_max = 0.00607504

    x_centers, x_edges = determine_1d_linear_bins(qx_min, qx_max, 5)
    y_centers, y_edges = determine_1d_linear_bins(qy_min, qy_max, 5)

    # verify
    gold_x_centers, gold_y_centers = get_gold_2d_linear_bins()
    gold_y_centers = gold_y_centers[::-1]

    np.allclose(x_centers, gold_x_centers, atol=5E-6)
    np.allclose(y_centers, gold_y_centers, atol=5E-6)

    # Check X
    assert abs(x_edges[1] - (-0.004694044)) < 1E-8
    assert abs(x_edges[2] - (-0.001814261)) < 1E-8
    # Check Y
    assert abs(y_edges[1] - (-0.002826)) < 1E-6
    assert abs(y_edges[2] - (-0.000601)) < 1E-6

    # Bin 2D
    # Get Q1D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test for no-weight binning
    binned_iq_2d = do_2d_no_weight_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                           x_edges, y_edges)

    # Verify Qx and Qy
    assert abs(x_centers[1] - (-0.003254)) < 1.E-6, 'Qx is not correct'
    assert abs(y_centers[1] - (-0.001713)) < 1.E-6, 'Qy is not correct'

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    assert abs(binned_iq_2d[0][1][1] - 67.) < 1E-6, 'I(Qx, Qy) is incorrect'
    assert abs(binned_iq_2d[1][1][1] - 4.725815626) < 1E-8, 'sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    # correct: 3.2999999999999996e-05
    assert abs(binned_iq_2d[2][1][1] - 3.31E-05) < 2E-7, 'dQx is incorrect'
    assert abs(binned_iq_2d[3][1][1] - 1.75E-05) < 2E-7, 'dQy is incorrect'

    # Test for weighted-binning
    binned_iq_2d = do_2d_weighted_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                          x_edges, y_edges)

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    # test value: 56.86602493293357
    assert abs(binned_iq_2d[0][1][1] - 56.8660) < 1E-4, 'Weighted-binned I(Qx, Qy) is incorrect'
    assert abs(binned_iq_2d[1][1][1] - 4.353773265) < 1E-8, 'Weighted-binned sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    assert abs(binned_iq_2d[2][1][1] - 3.30E-05) < 2E-7, 'dQx is incorrect'
    # correct: 1.71877860186208e-05
    assert abs(binned_iq_2d[3][1][1] - 1.75E-05) < 4E-7, 'dQy is incorrect'

    return


def get_gold_theta_bins():
    """Get theta bins from EXCEL

    Returns
    -------
    ndarray, ndarray
        bin edges, bin centers

    """
    theta_centers = np.ndarray([18, 54, 90, 126, 162, 198, 234, 270, 306, 342])
    theta_edges = np.ndarray([0, 36, 72, 108, 144, 180, 216, 252, 288, 324, 360])

    return theta_edges, theta_centers


def test_1d_annular_no_wt():
    """Test '1D_annular_no_sub_no_wt'

    Returns
    -------

    """
    theta_min = 0
    theta_max = 360.
    theta_step = 10.
    num_bins = int((theta_max - theta_min) / theta_step)

    theta_bin_centers, theta_bin_edges = determine_1d_linear_bins(theta_min, theta_max, num_bins)
    print(theta_bin_centers)
    print(theta_bin_edges)

    # Generate testing data: Get Q1D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Calculate theta array
    theta_array = np.arctan2(qy_array, qx_array) * 180. / np.pi
    # convert -0 to -180 to 180 to 360
    theta_array[np.where(theta_array < 0)] += 360.

    # calculate dQ from dQx and dQy
    dq_array = np.sqrt(dqx_array**2 + dqy_array**2)

    # binning
    binned_iq = do_1d_no_weight_binning(theta_array, dq_array, intensities, sigmas,
                                        theta_bin_centers, theta_bin_edges)

    # Check result
    print('Theta = 54 I[1]:  {} - {} = {}'.format(binned_iq.i[1], 63.66666667, binned_iq.i[1] - 63.66666667))
    print('Theta = 54 sI[1]: {} - {} = {}'.format(binned_iq.sigma[1], 3.257470048, binned_iq.sigma[1] - 3.257470048))

    # # Calculate azimuthal angle of each Q
    # angle = np.arctan2(qy, qx) * 180. / np.pi
    # angle[np.where(angle < 0)] += 360

    gold_theta_edges, gold_theta_centers = get_gold_theta_bins()

    assert np.allclose(theta_bin_centers, gold_theta_centers, rtol=1.e-5)
    assert np.allclose(theta_bin_edges, gold_theta_edges, rtol=1.e-5)

    return


def next3_test_1d_bin_log_wedge_no_wt():
    """Test '1D_bin_log_wedget_no_sub_no_wt

    Returns
    -------

    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    step_per_decade = 10  # 10 steps per decade

    bin_edges, bin_centers = determine_1d_log_bins(q_min, q_max, step_per_decade)

    # Bin wedge
    assert bin_edges
    assert bin_centers

    return
