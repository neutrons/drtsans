import numpy as np
from drtsans.mono.gpsans.prepare_sensitivity import prepare_sensitivity
import pytest

# This test implements issue #205 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
# DEV - Wenduo Zhou <petersonpf@ornl.gov> and Joe Osborn <osbornjd@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Lisa


def show_diff(test_data, gold_data):
    """

    Parameters
    ----------
    gold_data : ndarray
        gold data in 2D
    test_data : ndarray
        test data in 2D

    Returns
    -------
    None
    """
    if gold_data.shape != test_data.shape:
        # Data shape are different
        print('Gold data and test data are in different dimension: {} vs {}'
              ''.format(gold_data.shape, test_data.shape))

    else:
        print('      gold\ttest\tdiff')
        for m in range(gold_data.shape[0]):
            for n in range(gold_data.shape[1]):
                print('({}, {}):  {}\t{}\t{}'
                      .format(m, n, gold_data[m, n], test_data[m, n], abs(gold_data[m, n] - test_data[m, n])))

    return


def generate_test_data():
    """Generate test data from IS (Lisa)

    Returns
    -------
    np.ndarray, ndarray, ndarray, ndarray, ndarray, ndarray
        Matrix A, error of A, Matrix B, error of B, Matrix C, error of C
    """
    # Hard code Lisa's input in the Excel file
    l_matrix_a = np.array([[99, 400, 98],
                           [96, 95, np.nan],
                           [97, 98, 36]])

    l_matrix_b = np.array([[86, 400, 89],
                           [92, np.nan, 91],
                           [95, 97, 36]])

    l_matrix_c = np.array([[99, 400, 93],
                           [np.nan, 105, 94],
                           [95, 99, 36]])

    # Calculate uncertainties
    error_matrix_a = np.sqrt(l_matrix_a)
    error_matrix_b = np.sqrt(l_matrix_b)
    error_matrix_c = np.sqrt(l_matrix_c)

    return l_matrix_a, error_matrix_a, l_matrix_b, error_matrix_b, l_matrix_c, error_matrix_c


def get_weighted_averaged_matrix_abc():
    """Get the weighted average for 3 test files A, B, C with the uncertainties

    Returns
    -------
    ndarray, ndarray, ndarray, ndarray, ndarray, ndarray
    """
    m_a = np.array([
        [1.138985248, 4.601960598, 1.127480346],
        [1.104470543, 1.092965642, np.nan],
        [1.115975445, 1.127480346, 0.414176454]
    ])

    m_b = np.array([
        [1.0302349592, 4.7917905079, 1.066173388],
        [1.1021118168, np.nan, 1.0901323405],
        [1.1380502456, 1.1620091982, 0.4312611457]
    ])

    m_c = np.array([
        [1.137521253, 4.596045467, 1.068580571],
        [np.nan, 1.206461935, 1.080070685],
        [1.091560798, 1.137521253, 0.413644092]
    ])

    m_a_sigma = np.array([
        [0.122350147, 0.288793202, 0.121654003],
        [0.120254372, 0.119550801, np.nan],
        [0.120955425, 0.121654003, 0.070793758]
    ])

    m_b_sigma = np.array([
        [0.1180297776, 0.3029622398, 0.1203094805],
        [0.1225625632, np.nan, 0.1218144019],
        [0.1247904675, 0.1262624078, 0.0737887849],
    ])

    m_c_sigma = np.array([
        [0.122183097, 0.288354303, 0.11797513],
        [np.nan, 0.126304947, 0.118682824],
        [0.119387912, 0.122183097, 0.070700527]
    ])

    return m_a, m_a_sigma, m_b, m_b_sigma, m_c, m_c_sigma


def get_non_normalized_sensitivities():
    """Get the gold data for non-normalized pixel-wise sensitivities.

    Matrix D and sigmaD in Lisa's test document

    Returns
    -------

    """
    m_d = np.array([
        [1.1005340053, -np.inf, 1.0868702234],
        [1.1033135999, 1.1465981708, 1.0849705194],
        [1.1145185194, 1.1418683222, -np.inf],
    ])

    sigma_d = np.array([
        [0.0697462219, -np.inf, 0.0692534542],
        [0.0858371152, 0.0868247739, 0.0850069866],
        [0.0702338055, 0.0711963634, -np.inf],
    ])

    return m_d, sigma_d


def get_final_sensitivities():
    """Get the final sensitivities and sigma

    "Apply another weighted average and propate the error"

    Returns
    -------
    ndarray, ndarray

    """
    sen_matrix = np.array([
        [9.91E-01, -np.inf, 9.78E-01],
        [9.93E-01, 1.03E+00, 9.77E-01],
        [1.00E+00, 1.03E+00, -np.inf]
    ])

    sen_sigma_matrix = np.array([
        [1.63E-01, -np.inf, 1.62E-01],
        [1.70E-01, 1.76E-01, 1.67E-01],
        [1.65E-01, 1.69E-01, -np.inf]
    ])

    return sen_matrix, sen_sigma_matrix


def verify_result(sensitivity, sensitivity_error):
    """Verify the test result

    Parameters
    ----------
    sensitivity: ndarray
        sensitivity
    sensitivity_error: ndarray
        propagated sensitivity error

    Returns
    -------
    boolean
        result is same as gold value
    """
    assert sensitivity
    assert sensitivity_error


def normalize_by_monitor(flood_data, flood_data_error, monitor_counts):
    """Normalize the flood data field data by monitor

    Parameters
    ----------
    flood_data: ndarray
        flood data
    flood_data_error: ndarray
        flood data error
    monitor_counts: int/float
        monitor counts
    Returns
    -------
    ndarray, ndarray
        normalized flood data, normalized flood data error
    """
    return flood_data / monitor_counts, flood_data_error / monitor_counts


def calculate_weighted_average_with_error(normalized_data, normalized_error):
    """Calculated weighted average for normalized flood data and error

    Average = (sum I(m, n) / sigma^(m, n)) / sum 1 / sigma^2

    Parameters
    ----------
    normalized_data : ndarray
        normalized flood data
    normalized_error : ndarray
        normalized flood data's error

    Returns
    -------
    ndarray, ndarray, float, float
        data normalized by average, data's error normalized by average, Average, sigma(Average)

    """
    # Calculate weighted average
    # a = sum_{m, n} I(m, n) / sigma^2(m, n)
    weighted_sum = np.sum(normalized_data[~np.isnan(normalized_data)] /
                          normalized_error[~np.isnan(normalized_error)]**2)
    # b = sum 1 / sigma^2(m, n)
    weights_square = np.sum(1. / normalized_error[~np.isnan(normalized_error)]**2)
    # Avg = a / b
    weighted_average = weighted_sum / weights_square
    # sigma Avg = 1 / sqrt(b)
    weighted_average_error = 1. / np.sqrt(weights_square)
    print('[DEBUG] Average = {}, Error = {}'.format(weighted_average, weighted_average_error))

    # Normalize data by weighted-average
    avg_norm_data = normalized_data / weighted_average

    # Propagate uncertainties: sigma S(m, n) = I(m, n) / avg * [(error(m, n)/I(m, n))^2 + (sigma Avg/Avg)^2]^1/2
    avg_norm_error = normalized_data / weighted_average * np.sqrt((normalized_error / normalized_data)**2
                                                                  + (weighted_average_error / weighted_average)**2)

    return avg_norm_data, avg_norm_error, weighted_average, weighted_average_error


def process_bad_pixels(data, data_error, threshold_min, threshold_max):
    """Apply bad pixel threshold to each data set including error

    If any pixel with counts falling out of allowed threshold, i.e., out of range (min, max)
    they will be specified as bad pixels.
    For any bad pixel, the counts will then be set to '-inf'

    Parameters
    ----------
    data : ndarray
        normalized data
    data_error : ndarray
        normalized data error
    threshold_min: float
        minimum value of allowed value
    threshold_max: float
        maximum value of allowed value

    Returns
    -------
    ndarray, ndarray
        data with bad pixels set to INF,
        data error with bad pixels set to INF
    """
    data[(data < threshold_min) | (data > threshold_max)] = -np.inf
    data_error[(data < threshold_min) | (data > threshold_max)] = -np.inf

    print('[DEBUG] Process bad pixel:  Q threshold: {}, {}'.format(threshold_min, threshold_max))
    print(data)
    print('END')

    return data, data_error


def calculate_pixel_wise_sensitivity(data_a, data_a_error, data_b, data_b_error, data_c, data_c_error):
    """Calculate pixel-wise average of N files to create the new summed file for doing sensitivity correction

    D(m, n) = A_F(m, n) + B_F(m, n) + C_F(m, n) with average weight

    Calculate Pixel-wise Average of 3 files to create the new summed file for
    doing the sensitivity correction

    Returns
    -------
    nparray, nparray
        non-normalized sensitivities, non-normalized sensitivities error

    """
    # Create sensitivities and sigmas matrices
    sensitivities = np.zeros_like(data_a)
    sensitivities_error = np.zeros_like(data_a)

    # Calculate D'(i, j)    = sum_{k}^{A, B, C}M_k(i, j)/s_k^2(i, j)
    #           1/s^2(i, j) = sum_{k}^{A, B, C}1/s_k^2(i, j)
    for i in range(sensitivities.shape[0]):
        for j in range(sensitivities.shape[1]):
            d_ij_arr = np.array([data_a[i, j], data_b[i, j], data_c[i, j]])
            s_ij_arr = np.array([data_a_error[i, j], data_b_error[i, j], data_c_error[i, j]])

            if -np.inf in d_ij_arr:
                # Infinity case
                sensitivities[i, j] = -np.inf
                sensitivities_error[i, j] = -np.inf
            else:
                # Do weighted summation
                # remove NaN
                d_ij_arr = d_ij_arr[~np.isnan(d_ij_arr)]
                s_ij_arr = s_ij_arr[~np.isnan(s_ij_arr)]

                # sum
                s_ij = np.sum(1. / s_ij_arr**2)
                d_ij = np.sum(d_ij_arr / s_ij_arr**2) / s_ij
                s_ij = np.sqrt(s_ij)

                print('[DEBUG] ({}, {})  D" = {}, s = {}/{}  <-- {}, {}'
                      ''.format(i, j, d_ij, s_ij, s_ij**2, d_ij_arr, s_ij_arr))

                sensitivities[i, j] = d_ij
                sensitivities_error[i, j] = 1. / s_ij
    # END-FOR

    # Calculate  D, i.e., sensitivity for each pixel
    # D_avg(m, n) = sum_M^{a, b, c}{M(m, n) / sigma_M^2(m, n)} / sum_M^{A, B, C} 1 / sigma_M(m, n)**2
    #  sensitivities = (data_a / data_a_error**2 + data_b / data_b_error**2 + data_c / data_c_error**2) / \
    #                  (1 / data_a_error**2 + 1 / data_b_error**2 + 1 / data_c_error**2)

    # Debug output
    for i in range(sensitivities.shape[0]):
        row = ''
        for j in range(sensitivities.shape[1]):
            row += '{}\t'.format(sensitivities[i, j])
        print(row)
    # END-FOR

    # Debug output
    for i in range(sensitivities_error.shape[0]):
        row = ''
        for j in range(sensitivities_error.shape[1]):
            row += '{}\t'.format(sensitivities_error[i, j])
        print(row)
    # END-FOR

    return sensitivities, sensitivities_error


def normalize_sensitivities(d_matrix, sigma_d_matrix):
    """Do weighted average to pixel-wise sensitivities and propagate the error
    And then apply the average to sensitivity

    S_avg = sum_{m, n}{D(m, n) / sigma^2(m, n)} / sum_{m, n}{1 / sigma^2(m, n)}

    Parameters
    ----------
    d_matrix : ndarray
        pixel-wise sensitivities
    sigma_d_matrix : ndarray
        pixel-wise sensitivities error

    Returns
    -------
    ndarray, ndarray, float, float
        normalized pixel-wise sensitivities, normalized pixel-wise sensitivities error
        scalar sensitivity, error of scalar sensitivity
    """
    # Filter the matrix:
    dd_matrix = d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))]
    dd_sigma_matrix = sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))]
    print('Pure D matrix: {}'.format(dd_matrix))
    print('Pure sigma matrix: {}'.format(dd_sigma_matrix))

    # Calculate wighted-average of pixel-wise sensitivities: sum on (m, n)
    denomiator = np.sum(d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))] /
                        sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))]**2)
    nominator = np.sum(1 / sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))]**2)
    sens_avg = denomiator / nominator
    print('[DEBUG] S_avg = {} / {} = {}'.format(denomiator, nominator, sens_avg))

    # Normalize pixel-wise sensitivities
    sensitivities = d_matrix / sens_avg

    # Calculate scalar sensitivity's error
    # sigma_S_avg = sqrt(1 / sum_{m, n}(1 / sigma_D(m, n)^2))
    sigma_sens_avg = np.sqrt(1 / np.sum(1 / sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))]))

    # Propagate the sensitivities
    # sigma_sens(m, n) = D(m, n) / S_avg * [(sigma_D(m, n) / D(m, n))^2 + (sigma_S_avg / S_avg)^2]^{1/2}
    # D(m, n) are the non-normalized sensitivities
    sensitivities_error = d_matrix / sens_avg * np.sqrt((sigma_d_matrix / d_matrix)**2
                                                        + (sigma_sens_avg / sens_avg)**2)

    return sensitivities, sensitivities_error, sens_avg, sigma_sens_avg


def test_prepare_moving_det_sensitivity():
    """Test main algorithm to prepare sensitivity for instrument with moving detector

    Returns
    -------

    """
    # Set up the test data
    test_data_set = generate_test_data()
    monitor_a = 10
    monitor_b = 10
    monitor_c = 10
    threshold_min = 0.5
    threshold_max = 1.5

    # Normalize the flood field data by monitor: A, B and C
    matrix_a, sigma_a = test_data_set[0], test_data_set[1]
    matrix_a, sigma_a = normalize_by_monitor(matrix_a, sigma_a, monitor_a)

    matrix_b, sigma_b = test_data_set[2], test_data_set[3]
    matrix_b, sigma_b = normalize_by_monitor(matrix_b, sigma_b, monitor_b)

    matrix_c, sigma_c = test_data_set[4], test_data_set[5]
    matrix_c, sigma_c = normalize_by_monitor(matrix_c, sigma_c, monitor_c)

    # Find weighted average for each file and error
    matrix_a, sigma_a, avg_a, sigma_avg_a = calculate_weighted_average_with_error(normalized_data=matrix_a,
                                                                                  normalized_error=sigma_a)
    matrix_b, sigma_b, avg_b, sigma_avg_b = calculate_weighted_average_with_error(normalized_data=matrix_b,
                                                                                  normalized_error=sigma_b)
    matrix_c, sigma_c, avg_c, sigma_avg_c = calculate_weighted_average_with_error(normalized_data=matrix_c,
                                                                                  normalized_error=sigma_c)

    # check weighted sum
    # Matrix A
    assert abs(avg_a - 8.69) < 0.01, 'Average A is not correct'
    assert abs(sigma_avg_a - 0.329) < 0.001, 'Average Sigma_A is not correct'

    # Matrix C
    assert abs(avg_c - 8.7031340939) < 0.000001, 'Average C is not correct'
    assert abs(sigma_avg_c - 0.3298320424) < 0.0000001, 'Average Sigma_C is not correct'

    # Matrix B
    assert abs(avg_b - 8.347610342) < 0.00001, 'Average B is not correct: Test = {}, Gold = {}' \
                                               ''.format(avg_b, 7.3189956223)
    assert abs(sigma_avg_b - 0.323024967) < 0.000001, 'Average Sigma_B is not correct'

    # compare A, B and C to Lisa's gold data
    gold_data_set = get_weighted_averaged_matrix_abc()

    # A
    np.testing.assert_allclose(matrix_a, gold_data_set[0], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged A does not match gold data',
                               verbose=True)

    np.testing.assert_allclose(sigma_a, gold_data_set[1], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged sigma A does not match gold data',
                               verbose=True)

    # B
    np.testing.assert_allclose(matrix_b, gold_data_set[2], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged B does not match gold data',
                               verbose=True)
    np.testing.assert_allclose(sigma_b, gold_data_set[3], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged sigma B does not match gold data',
                               verbose=True)

    # C
    np.testing.assert_allclose(matrix_c, gold_data_set[4], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged C does not match gold data',
                               verbose=True)
    np.testing.assert_allclose(sigma_c, gold_data_set[5], rtol=1e-8, equal_nan=True,
                               err_msg='Weighted-averaged sigma C does not match gold data',
                               verbose=True)

    # Apply bad pixel threshold to each file
    matrix_a2, sigma_a2 = process_bad_pixels(matrix_a, sigma_a, threshold_min, threshold_max)
    matrix_b2, sigma_b2 = process_bad_pixels(matrix_b, sigma_b, threshold_min, threshold_max)
    matrix_c2, sigma_c2 = process_bad_pixels(matrix_c, sigma_c, threshold_min, threshold_max)

    #  there is no check for this step

    # Correct for beam stop and add N files for non-normalized sensitivities
    matrix_d, sigma_matrix_d = calculate_pixel_wise_sensitivity(matrix_a2, sigma_a2,
                                                                matrix_b2, sigma_b2,
                                                                matrix_c2, sigma_c2)

    # verify
    gold_d, gold_sigma_d = get_non_normalized_sensitivities()
    np.testing.assert_allclose(matrix_d, gold_d, rtol=1e-8, equal_nan=True,
                               err_msg='Summed files (matrix D) does not match to gold data',
                               verbose=True)
    np.testing.assert_allclose(sigma_matrix_d, gold_sigma_d, rtol=1e-8, equal_nan=True,
                               err_msg='Summed files (matrix sigma D) does not match to gold data',
                               verbose=True)

    # Normalize pixel-wise sensitivities by weighting-average
    sensitivities, sensitivities_error, avg_sens, avg_sigma_sens = normalize_sensitivities(matrix_d,
                                                                                           sigma_matrix_d)
    print('Avg sensitivities: {} - {} = {}'.format(avg_sens, 1.11E+00, abs(avg_sens - 1.11E+00)))

    gold_final_sen_matrix, gold_final_sigma_matrix = get_final_sensitivities()
    np.testing.assert_allclose(sensitivities, gold_final_sen_matrix, rtol=1e-2, equal_nan=True,
                               err_msg='Final sensitivities matrix not match', verbose=True)
    # FIXME : skip to wait for Lisa's corrected version
    # np.testing.assert_allclose(sensitivities_error, gold_final_sigma_matrix, rtol=1e-2, equal_nan=True,
    #                            err_msg='Final sensitivities error matrix not match', verbose=True)

    # Test high level method
    matrix_a, sigma_a = test_data_set[0], test_data_set[1]
    matrix_b, sigma_b = test_data_set[2], test_data_set[3]
    matrix_c, sigma_c = test_data_set[4], test_data_set[5]

    # convert input data to required format
    flood_matrix = np.ndarray(shape=(3, matrix_a.size), dtype=float)
    flood_error = np.ndarray(shape=(3, matrix_a.size), dtype=float)

    flood_matrix[0] = matrix_a.flatten()
    flood_matrix[1] = matrix_b.flatten()
    flood_matrix[2] = matrix_c.flatten()

    flood_error[0] = sigma_a.flatten()
    flood_error[1] = sigma_b.flatten()
    flood_error[2] = sigma_c.flatten()

    monitor_counts = np.array([monitor_a, monitor_b, monitor_c])

    test_sens_array, test_sens_sigma_array = prepare_sensitivity(flood_matrix, flood_error, monitor_counts,
                                                                 threshold_min, threshold_max)

    # # Test on Matrix A after normalized by weighted average and with bad pixels
    # np.testing.assert_allclose(matrix_a2.flatten(), test_sens_array[0], 1e-7)
    # np.testing.assert_allclose(sigma_a2.flatten(), test_sens_sigma_array[0], 1e-7)
    #
    # np.testing.assert_allclose(matrix_b2.flatten(), test_sens_array[1], 1e-7)
    # np.testing.assert_allclose(sigma_b2.flatten(), test_sens_sigma_array[1], 1e-7)
    #
    # np.testing.assert_allclose(matrix_c2.flatten(), test_sens_array[2], 1e-7)
    # np.testing.assert_allclose(sigma_c2.flatten(), test_sens_sigma_array[2], 1e-7)

    # Test raw sensitivities
    np.testing.assert_allclose(matrix_d.flatten(), test_sens_array, 1e-7)
    np.testing.assert_allclose(sigma_matrix_d.flatten(), test_sens_sigma_array, 1e-7)

    # np.testing.assert_allclose(sensitivities, test_sens_array, 1e-7)
    # np.testing.assert_allclose(sensitivities_error, test_sens_sigma_array, 1e-7)

    return


if __name__ == "__main__":
    pytest.main()
