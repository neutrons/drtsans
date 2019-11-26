import numpy as np

# This module supports testing data for issue #239, #245, #246 and #247.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/239
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/245
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/246
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/247
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data is from William's tests in eqsans_tof_q_binning_tests_R4.xlsx
# Intensities for a Generic 2D detector at 3 wave lengths
# Test EXCEL can be found at
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/5423db9b77dfd4911bf799c247530865/
# eqsans_tof_q_binning_tests_R4.xlsx

# The workspace is assued to have 3 wave lengths
intensities_matrix = np.array([[[93, 60, 89, 32, 97],
                                [43, 61, 82, 97, 55],
                                [78, 34, 50, 54, 67],
                                [98, 88, 37, 92, 97],
                                [72, 97, 100, 71, 39]],  # 3.0 A
                               [[76, 39, 51, 70, 61],
                                [64, 54, 78, 35, 30],
                                [67, 98, 100, 56, 79],
                                [97, 35, 41, 90, 45],
                                [30, 41, 68, 34, 51]],  # 3.1 A
                               [[78, 36, 46, 75, 91],
                                [64, 56, 92, 73, 60],
                                [74, 72, 69, 84, 87],
                                [36, 78, 40, 68, 72],
                                [59, 40, 39, 34, 85]]  # 3.2A
                               ])

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

# Scalar dQ matrix copied from revision 3's 1D_bin_linear_no_sub_no_wt and 1D_bin_log_no_sub_no_wt
scalar_dq_matrix = np.array([
    # 3.0 A
    [[0.011912, 0.011912, 0.011912, 0.011912, 0.011912],
     [0.011912, 0.011912, 0.011912, 0.011912, 0.011912],
     [0.011912, 0.011912, 0.011912, 0.011912, 0.011912],
     [0.011912, 0.011912, 0.011912, 0.011912, 0.011912],
     [0.011912, 0.011912, 0.011912, 0.011912, 0.011912]],
    # 3.1 A
    [[0.011527, 0.011527, 0.011527, 0.011527, 0.011527],
     [0.011527, 0.011527, 0.011527, 0.011527, 0.011527],
     [0.011527, 0.011527, 0.011527, 0.011527, 0.011527],
     [0.011527, 0.011527, 0.011527, 0.011527, 0.011527],
     [0.011527, 0.011527, 0.011527, 0.011527, 0.011527]],
    # 3.2 A
    [[0.011167, 0.011167, 0.011167, 0.011167, 0.011167],
     [0.011167, 0.011167, 0.011167, 0.011167, 0.011167],
     [0.011167, 0.011167, 0.011167, 0.011167, 0.011167],
     [0.011167, 0.011167, 0.011167, 0.011167, 0.011167],
     [0.011167, 0.011167, 0.011167, 0.011167, 0.011167]]
])


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
    # Check input: dimension must be either 1 or 2
    if q_dimension not in [1, 2]:
        raise RuntimeError('Q-dimension must be 1 or 2')

    if q_dimension == 1:
        # Calculate scalar Q
        scalar_q_matrix = np.sqrt(qx_matrix ** 2 + qy_matrix ** 2)
        # Scalar dQ is defined
    else:
        # No-op
        scalar_q_matrix = None

    # Define a None returning object
    returns = None
    if drt_standard:
        # Convert the Q, dQ, I, sigma I in 2D matrix to 1D arrays to match drtsans binning methods' requirements
        intensity_array = intensities_matrix.flatten()
        sigma_array = uncertainties_matrix.flatten()
        if q_dimension == 1:
            # Q1D: scalar Q and scalar dQ from matrix to 1D array
            scalar_q_array = scalar_q_matrix.flatten()
            scalar_dq_array = scalar_dq_matrix.flatten()
            returns = intensity_array, sigma_array, scalar_q_array, scalar_dq_array
        else:
            # Q2D: Qx, Qy, dQx and dQy from matrix to 1D array
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


def get_gold_wedge_angles():
    """Get wedge angles from Excel for test

    Returns
    -------
    ndarray
        Wedge angles for all the pixels

    """
    wedge_angle_array = np.array([
        141.026949, 123.2554061, 94.31432613, 63.21170404, 42.66018375,
        155.9524823, 139.9324658, 97.78839725, 47.53052642, 26.94261266,
        175.230287, 171.0616985, 126.183895, 11.5457636, 5.429158913,
        195.6073192, 207.7689884, 257.6752727, -34.3685385, -17.65012424,
        212.7055538, 230.4368877, 264.5704513, -57.539898, -36.18192273,
    ])

    return wedge_angle_array


def get_gold_theta_bins():
    """Get theta bins from EXCEL

    Returns
    -------
    ndarray, ndarray
        bin edges, bin centers

    """
    theta_centers = np.array([18, 54, 90, 126, 162, 198, 234, 270, 306, 342])
    theta_edges = np.array([0, 36, 72, 108, 144, 180, 216, 252, 288, 324, 360])

    return theta_edges, theta_centers


def get_gold_azimuthal_values():
    """Get the azimuthal values from EXCEL

    Returns
    -------

    """
    gold_theta_array = np.array([141.026949, 123.2554061, 94.31432613, 63.21170404, 42.66018375,
                                 155.9524823, 139.9324658, 97.78839725, 47.53052642, 26.94261266,
                                 175.230287, 171.0616985, 126.183895, 11.5457636, 5.429158913,
                                 195.6073192, 207.7689884, 257.6752727, 325.6314615, 342.3498758,
                                 212.7055538, 230.4368877, 264.5704513, 302.460102, 323.8180773])

    return gold_theta_array


