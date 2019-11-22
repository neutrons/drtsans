import numpy as np
from drtsans.dataobjects import IQazimuthal, IQmod
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/iq.py
from drtsans.iq import _determine_1d_linear_bins, _determine_1d_log_bins, _do_1d_no_weight_binning,\
    BinningMethod, _do_2d_weighted_binning, _do_2d_no_weight_binning,\
    bin_annular_into_q1d, bin_wedge_into_q1d, BinningParams, bin_intensity_into_q1d, bin_iq_into_linear_q2d
import pytest

# This test implements issue #169 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/tree/169_bin_q1d
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data is from William's tests in eqsans_tof_q_binning_tests_R3.xlsx
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


def test_1d_bin_linear_no_wt():
    """Test case: '1D_bin_linear_no_sub_wt'

    Test methods for 1D linear no-weight binning

    Returns
    -------
    None

    """
    # From Tab '1D_bin_linear_no_sub_wt'
    q_min = 0.000
    q_max = 0.010
    num_bins = 10

    # Verify bin edges and bin center
    bin_centers, bin_edges = _determine_1d_linear_bins(q_min, q_max, num_bins)
    gold_edges, gold_centers = get_gold_1d_linear_bins()

    np.testing.assert_allclose(bin_edges, gold_edges, rtol=1.E-12)
    np.testing.assert_allclose(bin_centers, gold_centers, rtol=1.E-12)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = _do_1d_no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                         bin_centers, bin_edges)

    # Calculate and verify
    # I(0.0035) = 68.92857:    drtsans: 68.92857142857143
    # verify Q[3]
    assert abs(binned_iq.mod_q[3] - 0.0035) < 1E-6, 'Q[3] {} shall be {} +/- 1e-6' \
                                                    ''.format(binned_iq.delta_mod_q[3], 0.0035)
    # verify I[3]
    assert abs(binned_iq.intensity[3] - 68.92857) < 1E-5, 'Intensity[3] shall be 68.92857 but not {}' \
                                                          ''.format(binned_iq.intensity[3])
    # verify sigmaI[3] = 2.218889:
    assert abs(binned_iq.error[3] - 2.218889) < 1E-6, 'error'
    # verify sigma_Q[3] = 1.154E-02
    assert binned_iq.delta_mod_q[3] == pytest.approx(1.154e-02, abs=2.E-5), \
        'Linear binning: Q resolution {} does not match expected {}'.format(binned_iq.delta_mod_q[3], 1.135E-02)

    # Test high level method
    test_iq = IQmod(intensities, sigmas, scalar_q_array, scalar_dq_array, None)
    binning = BinningParams(q_min, q_max, num_bins)
    binned_iq2 = bin_intensity_into_q1d(test_iq, binning, True, BinningMethod.NOWEIGHT)

    # verify
    np.testing.assert_allclose(binned_iq2.intensity, binned_iq.intensity, rtol=1e-8,
                               equal_nan=True, err_msg='High level method cannot have same result from low levels',
                               verbose=True)

    return


def test_1d_bin_log_no_wt():
    """Test '1D_bin_log_no_sub_no_wt'

    Test binning methods for 1D no-weight binning with log bins

    Returns
    -------

    """
    # Define Q range from tab '1D_bin_log_no_sub_no_wt' in r4
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    num_steps_per_10 = 10  # 10 steps per decade

    # Verify bin edges and bin center
    bin_centers, bin_edges = _determine_1d_log_bins(q_min, q_max, num_steps_per_10)
    gold_edges, gold_centers = get_gold_1d_log_bins()

    np.testing.assert_allclose(bin_edges, gold_edges, rtol=5.E-4)
    np.testing.assert_allclose(bin_centers, gold_centers, rtol=5.E-4)

    # Get Q1D data
    intensities, sigmas, scalar_q_array, scalar_dq_array = generate_test_data(1, True)

    # Binned I(Q) no-weight
    binned_iq = _do_1d_no_weight_binning(scalar_q_array, scalar_dq_array, intensities, sigmas,
                                         bin_centers, bin_edges)

    # Verify: 2 I(Q) in bin: Q(3, 2, 3.1), Q(3, 2, 3.2)
    # I(0.0022) = 70.00000
    assert binned_iq.intensity[3] == pytest.approx(70.00000, abs=1.E-12), 'intensity'
    # dI(0.0022) = 5.9160797831
    assert binned_iq.error[3] == pytest.approx(5.9160797831, abs=1.E-12), 'error'
    # sigma_Q(0.0022) = 1.135E-02
    assert binned_iq.delta_mod_q[3] == pytest.approx(1.135E-02, abs=2.E-5),\
        'Log binning: Q resolution {} does not match expected {}'.format(binned_iq.delta_mod_q[3], 1.135E-02)

    # Test the high level method
    test_iq = IQmod(intensities, sigmas, scalar_q_array, scalar_dq_array)
    binning = BinningParams(q_min, q_max, num_steps_per_10)
    binned_iq = bin_intensity_into_q1d(test_iq, binning, False, BinningMethod.NOWEIGHT)
    assert binned_iq.intensity[3] == pytest.approx(70.00000, abs=1.E-12), 'intensity'


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

    x_centers, x_edges = _determine_1d_linear_bins(qx_min, qx_max, 5)
    y_centers, y_edges = _determine_1d_linear_bins(qy_min, qy_max, 5)

    # verify
    gold_x_centers, gold_y_centers = get_gold_2d_linear_bins()
    gold_y_centers = gold_y_centers[::-1]

    np.testing.assert_allclose(x_centers, gold_x_centers, atol=5E-6)
    np.testing.assert_allclose(y_centers, gold_y_centers, atol=5E-6)

    # Check X
    assert x_edges[1] == pytest.approx(-0.004694044, abs=1E-8)
    assert x_edges[2] == pytest.approx(-0.001814261, abs=1E-8)
    # Check Y
    assert y_edges[1] == pytest.approx(-0.002826, abs=1E-6)
    assert y_edges[2] == pytest.approx(-0.000601, abs=1E-6)

    # Bin 2D No-weight
    # Get Q1D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test for no-weight binning
    binned_iq_2d = _do_2d_no_weight_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                            x_edges, y_edges)

    # Verify Qx and Qy
    assert x_centers[1] == pytest.approx(-0.003254, abs=1.E-6), 'Qx is not correct'
    assert y_centers[1] == pytest.approx(-0.001713, abs=1.E-6), 'Qy is not correct'

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    assert binned_iq_2d[0][1][1] == pytest.approx(67., abs=1E-6), 'I(Qx, Qy) is incorrect'
    assert binned_iq_2d[1][1][1] == pytest.approx(4.725815626, abs=1E-8), 'sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    # correct: 3.2999999999999996e-05
    assert binned_iq_2d[2][1][1] == pytest.approx(3.31E-05, abs=2E-7), 'dQx is incorrect'
    assert binned_iq_2d[3][1][1] == pytest.approx(1.75E-05, abs=2E-7), 'dQy is incorrect'

    # Bin 2D Weighted
    # Test for weighted-binning
    binned_iq_2d = _do_2d_weighted_binning(qx_array, dqx_array, qy_array, dqy_array, intensities, sigmas,
                                           x_edges, y_edges)

    # verify I(-0.003254,-0.001713) and sigma(-0.003254,-0.001713)
    # test value: 56.86602493293357
    assert binned_iq_2d[0][1][1] == pytest.approx(56.8660, abs=1E-4), 'Weighted-binned I(Qx, Qy) is incorrect'
    assert binned_iq_2d[1][1][1] == pytest.approx(4.353773265, abs=1E-8), \
        'Weighted-binned sigma I(Qx, Qy) is incorrect'

    # verify dQx and dQy
    assert binned_iq_2d[2][1][1] == pytest.approx(3.30E-05, abs=2E-7), 'dQx is incorrect'
    # correct: 1.71877860186208e-05
    assert binned_iq_2d[3][1][1] == pytest.approx(1.75E-05, abs=4E-7), 'dQy is incorrect'

    # Test API for high level method
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)
    qx_bin_params = BinningParams(qx_min, qx_max, 5)
    qy_bin_params = BinningParams(qy_min, qy_max, 5)
    binned_iq_2d_wt = bin_iq_into_linear_q2d(test_i_q, qx_bin_params, qy_bin_params, BinningMethod.WEIGHTED)
    # verify
    np.testing.assert_allclose(binned_iq_2d[0], binned_iq_2d_wt.intensity, atol=1E-10)

    return


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


# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/246
def test_1d_annular_no_wt():
    """Test '1D_annular_no_sub_no_wt'

    Returns
    -------

    """
    theta_min = 0
    theta_max = 360.
    num_bins = 10

    q_min = 0.003
    q_max = 0.006

    theta_bin_centers, theta_bin_edges = _determine_1d_linear_bins(theta_min, theta_max, num_bins)

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

    # binning
    binned_iq = _do_1d_no_weight_binning(theta_array[allowed_q_index],
                                         dq_array[allowed_q_index],
                                         intensities[allowed_q_index],
                                         sigmas[allowed_q_index],
                                         theta_bin_centers, theta_bin_edges)

    # Check bins
    gold_theta_edges, gold_theta_centers = get_gold_theta_bins()
    np.testing.assert_allclose(theta_bin_centers, gold_theta_centers, rtol=1.e-5)
    np.testing.assert_allclose(theta_bin_edges, gold_theta_edges, rtol=1.e-5)

    # Check theta (azimuthal angle)
    # print(theta_array)
    # assert abs(theta_array[0] - 141.026949) < 5E-2, 'Azimuthal angle check'
    gold_theta_array = get_gold_azimuthal_values()
    num_test_data = gold_theta_array.shape[0]

    np.testing.assert_allclose(theta_array[:num_test_data], gold_theta_array, rtol=6.e-4, err_msg='Azimuthal vectors')

    assert binned_iq.intensity[1] == pytest.approx(63.66666667, abs=1E-8), 'Binned intensity is wrong'
    assert binned_iq.error[1] == pytest.approx(3.257470048, abs=1E-8), 'Binned sigma I is wrong'
    # 4.70549611605334e-05 calculated vs 4.717e-05
    assert binned_iq.delta_mod_q[1] == pytest.approx(4.717E-05, abs=1.5E-7), 'Binned Q resolution is wrong'

    # Test the high level method
    # Define input data
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)

    # Bin
    theta_binning = BinningParams(theta_min, theta_max, num_bins)
    binned_iq = bin_annular_into_q1d(test_i_q, theta_binning, q_min, q_max, BinningMethod.NOWEIGHT)
    # verify
    assert binned_iq.intensity[1] == pytest.approx(63.66666667, abs=1E-8), 'Binned intensity is wrong'
    assert binned_iq.error[1] == pytest.approx(3.257470048, abs=1E-8), 'Binned sigma I is wrong'
    # 4.70549611605334e-05 calculated vs 4.717e-05
    assert binned_iq.delta_mod_q[1] == pytest.approx(4.717E-05, abs=1.5E-7), 'Binned Q resolution is wrong'


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


def test_1d_bin_log_wedge_no_wt():
    """Test '1D_bin_log_wedget_no_sub_no_wt
    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    step_per_decade = 10  # 10 steps per decade

    min_wedge_angle = -45.
    max_wedge_angle = 45

    # Bin wedge
    bin_centers, bin_edges = _determine_1d_log_bins(q_min, q_max, step_per_decade)

    # Get data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)
    # calculate Q and dQ
    scalar_q_array = np.sqrt(qx_array**2 + qy_array**2)
    scalar_dq_array = np.sqrt(dqx_array**2 + dqy_array**2)

    # Calculate wedge angles for each I(Qx, Qy)
    # calculate azimuthal angles from -180 to 180 degrees
    azimuthal_array = np.arctan2(qy_array, qx_array) * 180. / np.pi
    # correct azimuthal angles to -90 to 270 degrees
    azimuthal_array[azimuthal_array < -90.] += 360.

    # Define the filter (mask/ROI) for pixels falling into preferred wedge
    wedge_indexes = (azimuthal_array > min_wedge_angle) & (azimuthal_array < max_wedge_angle)

    # Binning
    binned_iq = _do_1d_no_weight_binning(q_array=scalar_q_array[wedge_indexes],
                                         dq_array=scalar_dq_array[wedge_indexes],
                                         iq_array=intensities[wedge_indexes],
                                         sigmaq_array=sigmas[wedge_indexes],
                                         bin_centers=bin_centers,
                                         bin_edges=bin_edges)

    # Verification
    # Bin centers and boundaries
    gold_edges, gold_centers = get_gold_1d_log_bins()
    np.testing.assert_allclose(bin_edges, gold_edges, rtol=5.E-4)
    np.testing.assert_allclose(bin_centers, gold_centers, rtol=5.E-4)

    # Azimuthal angles
    gold_angles = get_gold_wedge_angles()
    np.testing.assert_allclose(azimuthal_array[:gold_angles.shape[0]], gold_angles,
                               rtol=1e-3, err_msg='Wedge angles do no match to gold data')

    # Number of pixels in the wedge
    assert intensities[wedge_indexes].shape[0] == 7 * 3, 'Number of I(Q) in wedge area is incorrect'

    # Binned I(Q) and others - these don't match the values in the test below
    print('Q = 0.005623 [7]:  (Q, I, sigmaI, dQ)'
          '\nTest    : {}\t{}\t{}\t{}'
          '\nExpected: {}\t{}\t{}\t{}'.format(binned_iq.mod_q[7], binned_iq.intensity[7],
                                              binned_iq.error[7], binned_iq.delta_mod_q[7],
                                              0.005623, 65.77777778, 2.703450013, 5.798E-05))

    # NOTE: using wzz's "correct" value till William's value
    """
    Q = 0.005623 [7]:  (Q, I, sigmaI, dQ)
    Test    : 0.005623413251903491, 67.7, 2.6019223662515376, 5.8479318657713166e-05
    """
    assert binned_iq.intensity[7] == pytest.approx(67.7, abs=1E-10)
    assert binned_iq.error[7] == pytest.approx(2.601922, abs=1E-5)
    assert binned_iq.delta_mod_q[7] == pytest.approx(5.84793186e-05, abs=1E-10)

    # Test high level method
    # Define input data
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)
    binning = BinningParams(q_min, q_max, step_per_decade)

    binned_iq2 = bin_wedge_into_q1d(test_i_q, min_wedge_angle, max_wedge_angle,
                                    binning, linear_binning=False,
                                    method=BinningMethod.NOWEIGHT)

    # verify calculated I, sigma and dQ
    assert binned_iq.intensity[7] == pytest.approx(binned_iq2.intensity[7], abs=1E-12)
    assert binned_iq.error[7] == pytest.approx(binned_iq2.error[7], abs=1E-12)
    assert binned_iq.delta_mod_q[7] == pytest.approx(binned_iq2.delta_mod_q[7], abs=1E-12)


if __name__ == '__main__':
    pytest.main([__file__])
