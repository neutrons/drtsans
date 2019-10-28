import numpy as np
from drtsans.iq import IofQ

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
    # blabla

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
    q_min = 0.
    q_max = 0.010
    num_bins = 10

    # Verify bin edges and bin center
    bin_centers, bin_edges = determine_linear_bins(q_min, q_max, num_bins)
    gold_edges, gold_centers = get_gold_1d_linear_bins()

    assert np.allclose(bin_edges, gold_edges, 1.E-12)
    assert np.allclose(bin_centers, gold_centers, 1.E-12)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                  bin_centers, bin_edges)

    # Calculate and verify
    # I(0.0035) =		68.92857
    assert abs(binned_iq.i[3] - 68.92857) < 1.E-12, 'I wrong'
    # di(0.0035)		2.218889
    assert abs(binned_iq.sigma[3] - 2.218889) < 1.E-12, 'sigma I wrong'
    # sigma_Q(0.0035) = 		3.722E-05
    assert abs(binned_iq.dq[3] - 3.722E-05) < 1.E-12, 'Q resolution wrong'

    return


def next_test_1d_bin_log_no_wt():
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
    bin_edges, bin_centers = determine_1d_log_bins(q_min, q_max, num_steps_per_10)
    gold_edges, gold_centers = get_gold_1d_log_bins()

    assert np.allclose(bin_edges, gold_edges, 1.E-12)
    assert np.allclose(bin_centers, gold_centers, 1.E-12)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                  bin_centers, bin_edges)

    # Verify: 2 I(Q) in bin: Q(3, 2, 3.1), Q(3, 2, 3.2)
    # I(0.0022) = 70.00000
    assert abs(binned_iq.i[4] - 70.00000) < 1.E-12, 'I wrong'
    # dI(0.0022) = 5.9160797831
    assert abs(binned_iq.sigma[4] - 5.9160797831) < 1.E-12, 'sigma I wrong'
    # sigma_Q(0.0022) = 2.529E-05
    assert abs(binned_iq.dq[4] - 2.529E-05) < 1.E-12, 'Q resolution wrong'

    return


def next_test_1d_bin_log_wedge_no_wt():
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


def next3_test_1d_annular_no_wt():
    """Test '1D_annular_no_sub_no_wt'

    Returns
    -------

    """

    return


def next2test_2d_bin_no_wt():
    """Test '2D_bin_no_sub_no_wt'

    Returns
    -------

    """

    return


def determine_2d_bins():
    """

    Returns
    -------

    """
    print('Hello World')

    min_q = -0.006134
    min_q_near = - 0.003254

    qx_min_edge = min_q - 0.5 * (min_q_near - min_q)
    qx_max_edge = 0.0068250907

    qx_step = (qx_max_edge - qx_min_edge) / 5

    print(qx_min_edge)
    print(qx_step)
    print(qx_step - 0.002880)

    return


# TODO FIXME - Move this method back to drtsans.iq
def determine_linear_bins(q_min, q_max, bins):
    """Determine linear bin edges and centers

    Parameters
    ----------
    q_min : float
        Q min of bin edge
    q_max : float
        Q max of bin edge
    bins : integer
        number of bins

    Returns
    -------
    ndarray, ndarray
        bin centers, bin edges
    """
    delta_q = (q_max - q_min) / bins
    bin_edges = np.arange(bins + 1).astype('float') * delta_q + q_min
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

    return bin_centers, bin_edges


# TODO FIXME - move this method back to drtsans.iq with different name to Lisa's
def determine_1d_log_bins(q_min, q_max, step_per_decade):
    """

    Parameters
    ----------
    q_min
    q_max
    step_per_decade: float
        step per decade (ex. 0.1 to 1.0 is one decade); denoted as 'j' in document
    Returns
    -------

    """
    # Calculate step and align q_min to q0, a decade (power of 10) nearest to q_min but less than q_min
    # 20191016 IS: "0.2% error is allowed.  This formula ensure that the number of steps per decade is respected"
    delta = np.power(10., 1. / step_per_decade)
    q0 = np.power(delta, np.floor(step_per_decade * np.log10(q_min)))
    print('[DEBUG OUTPUT: q_min = {}, q0 = {}'.format(q_min, q0))

    # Determine number of bins
    num_bins = 1 + int(np.ceil(step_per_decade * np.log(q_max / q0) / np.log(10)))
    print('[DEBUG OUTPUT: number of bins = {}'.format(num_bins))

    # Calculate bin centers
    bin_centers = np.arange(num_bins).astype('float')
    bin_centers = q0 * np.power(delta, bin_centers)

    # Calculate bin boundaries
    delta_q_array = 2. * (delta - 1) / (delta + 1) * bin_centers
    bin_edges = np.zeros((num_bins + 1,), dtype='float')
    bin_edges[1:] = bin_centers[:] + 0.5 * delta_q_array[:]
    bin_edges[0] = bin_centers[0] - 0.5 * delta_q_array[0]

    # # Big debug
    # print('[DEBUG OUTPUT] Edge from {}'.format(bin_edges[0]))
    # for i in range(99):
    #     from_left = bin_centers[i] + 0.5 * delta_q_array[i]
    #     from_right = bin_centers[i + 1] - 0.5 * delta_q_array[i + 1]
    #     diff = from_right - from_left
    #     diff2 = from_right - bin_edges[i+1]
    #     print('[DEBUG OUTPUT] From left = {}, From right = {}, Difference = {}, Calculated = {}  Diff = {}'
    #           ''.format(from_left, from_right, diff, bin_edges[i+1], diff2))
    # # END-FOR

    return bin_centers, bin_edges


def weighted_binning(q_array, dq_array, iq_array, sigma_iq_array, bin_centers, bin_edges):
    """ Bin I(Q) by given bin edges and do weighted binning

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

    Parameters
    ----------
    q_array: ndarray
        scaler momentum transfer Q
    dq_array: ndarray
        scaler momentum transfer (Q) resolution
    iq_array: ndarray
        I(Q)
    sigma_iq_array: ndarray
        sigma I(Q)
    bin_centers: numpy.ndarray
        bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
    bin_edges: numpy.ndarray
        bin edges
    Returns
    -------
    IofQ
        named tuple for Q, dQ, binned I(Q), binned sigma_I(Q)
    """
    # check input
    assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

    # Flatten input data to 1D
    # q_array = IofQCalculator.flatten(q_array)
    # dq_array = IofQCalculator.flatten(dq_array)
    # iq_array = IofQCalculator.flatten(iq_array)
    # sigma_iq_array = IofQCalculator.flatten(sigma_iq_array)

    # calculate 1/sigma^2 for multiple uses
    invert_sigma2_array = 1. / (sigma_iq_array ** 2)

    # Counts per bin: I_{k, raw} = \sum \frac{I(i, j)}{(\sigma I(i, j))^2}
    i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array * invert_sigma2_array)

    # Weight per bin: w_k = \sum \frac{1}{\sqrt{I(i, j)^2}
    w_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=invert_sigma2_array)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{w_k}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_array / w_array
    sigma_final_array = 1 / np.sqrt(w_array)

    # Calculate Q resolution of binned
    # FIXME - waiting for Lisa's equations for binned q resolution
    # FIXME - this is an incorrect solution temporarily for workflow
    binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
    bin_q_resolution = binned_dq / i_raw_array

    # Get the final result
    binned_iq = IofQ(bin_centers, bin_q_resolution, i_final_array, sigma_final_array)

    return binned_iq


def no_weight_binning(q_array, dq_array, iq_array, sigmaq_array, bin_centers, bin_edges):
    """ Bin I(Q) by given bin edges and do no-weight binning
    This method implements equation 11.34, 11.35 and 11.36 in master document.
    Parameters
    ----------
    q_array: ndarray
        scaler momentum transfer Q
    dq_array: ndarray
        scaler momentum transfer (Q) resolution
    iq_array: ndarray
        I(Q)
    sigmaq_array: ndarray
        sigma I(Q)
    bin_centers: numpy.ndarray
        bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
    bin_edges: numpy.ndarray
        bin edges
    Returns
    -------
    IofQ
        named tuple for Q, dQ, binned I(Q), binned sigma_I(Q)
    """
    # check input
    assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

    # Number of I(q) in each target Q bin
    num_pt_array, bin_x = np.histogram(q_array, bins=bin_edges)

    # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
    i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array)
    # Square of summed uncertainties for each bin
    sigma_sqr_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=sigmaq_array ** 2)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{Nk}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_array / num_pt_array
    sigma_final_array = np.sqrt(sigma_sqr_array) / num_pt_array

    # Calculate Q resolution of binned
    # FIXME - waiting for Lisa's equations for binned q resolution
    # FIXME - this is an incorrect solution temporarily for workflow
    binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
    bin_q_resolution = binned_dq / num_pt_array

    # Get the final result
    binned_iq = IofQ(bin_centers, bin_q_resolution, i_final_array, sigma_final_array)

    return binned_iq
