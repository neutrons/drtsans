import numpy as np
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
# from mantid.simpleapi import LoadEmptyInstrument, AddSampleLog

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
                print('({}, {}):  {}\t{}\t'
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
        [1.175024613, 5.465230759, 1.216013844],
        [1.257003075, np.nan, 1.243339998],
        [1.297992305, 1.325318459, 0.491870768]
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
        [0.136908362, 0.364490799, 0.139623949],
        [0.142310717, np.nan, 0.141418228],
        [0.144970268, 0.146728914, 0.084804125]
    ])

    m_c_sigma = np.array([
        [0.122183097, 0.288354303, 0.11797513],
        [np.nan, 0.126304947, 0.118682824],
        [0.119387912, 0.122183097, 0.070700527]
    ])

    return m_a, m_a_sigma, m_b, m_b_sigma, m_c, m_c_sigma


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
    # Calculate  D, i.e., sensitivity for each pixel
    # D_avg(m, n) = sum_M^{a, b, c}{M(m, n) / sigma_M^2(m, n)} / sum_M^{A, B, C} 1 / sigma_M(m, n)**2
    sensitivities = (data_a / data_a_error**2 + data_b / data_b_error**2 + data_c / data_c_error**2) / \
                    (1 / data_a_error**2 + 1 / data_b_error**2 + 1 / data_c_error**2)

    # Propagating error to D: delta Sen(m, n) = 1 / sum_{M}^{A, B, C}{ 1 / sigma_M^2(m, n)}
    sensitivities_error = 1 / (1 / data_a_error**2 + 1 / data_b_error*2 + 1 / data_c_error**2)

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
    # Calculate wighted-average of pixel-wise sensitivities: sum on (m, n)
    denomiator = np.sum(d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))] /
                        sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))])
    nominator = np.sum(1 / sigma_d_matrix[~(np.isinf(d_matrix) | np.isnan(d_matrix))])
    sens_avg = denomiator / nominator

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

    # Normalize by monitor: A, B and C
    matrix_a, sigma_a = test_data_set[0], test_data_set[1]
    matrix_a, sigma_a = normalize_by_monitor(matrix_a, sigma_a, monitor_a)

    matrix_b, sigma_b = test_data_set[2], test_data_set[3]
    matrix_b, sigma_b = normalize_by_monitor(matrix_b, sigma_b, monitor_b)

    matrix_c, sigma_c = test_data_set[4], test_data_set[5]
    matrix_c, sigma_c = normalize_by_monitor(matrix_c, sigma_c, monitor_c)

    # Find weighted average and error
    matrix_a, sigma_a, avg_a, sigma_avg_a = calculate_weighted_average_with_error(normalized_data=matrix_a,
                                                                                  normalized_error=sigma_a)
    matrix_b, sigma_b, avg_b, sigma_avg_b = calculate_weighted_average_with_error(normalized_data=matrix_b,
                                                                                  normalized_error=sigma_b)
    matrix_c, sigma_c, avg_c, sigma_avg_c = calculate_weighted_average_with_error(normalized_data=matrix_c,
                                                                                  normalized_error=sigma_c)

    assert abs(avg_a - 8.69) < 0.01, 'Average A is not correct'
    assert abs(sigma_avg_a - 0.329) < 0.001, 'Average Sigma_A is not correct'

    # compare A, B and C to Lisa's gold data
    gold_data_set = get_weighted_averaged_matrix_abc()

    # A
    if not np.allclose(matrix_a[~np.isnan(matrix_a)],
                       gold_data_set[0][~np.isnan(gold_data_set[0])], 1E-8):
        show_diff(matrix_a, gold_data_set[0])
        assert False, 'Weighted-averaged A does not match gold data'
    if not np.allclose(sigma_a[~np.isnan(sigma_a)],
                       gold_data_set[1][~np.isnan(gold_data_set[1])], 1E-8):
        show_diff(sigma_a, gold_data_set[1])
        assert False, 'Weighted-averaged sigma A does not match gold data'

    # B
    if not np.allclose(matrix_b, gold_data_set[2], 1E-8, True):
        print('B are different')
        show_diff(matrix_b, gold_data_set[2])
        # assert False, 'Weighted-averaged B does not match gold data'
    if not np.allclose(sigma_b, gold_data_set[3], 1E-8, True):
        print('Sigma B are different')
        show_diff(sigma_b, gold_data_set[3])
        # assert False, 'Weighted-averaged sigma B does not match gold data'

    # A
    if not np.allclose(matrix_c, gold_data_set[4], 1E-8, True):
        print('C are different')
        show_diff(matrix_c, gold_data_set[4])
        # assert False, 'Weighted-averaged A does not match gold data'
    if not np.allclose(sigma_c, gold_data_set[5], 1E-8, True):
        print('sigma C are different')
        show_diff(sigma_c, gold_data_set[5])
        assert False, 'Weighted-averaged sigma C does not match gold data'

    # Apply bad pixel threshold to the data
    matrix_a, sigma_a = process_bad_pixels(matrix_a, sigma_a, threshold_min, threshold_max)
    matrix_b, sigma_b = process_bad_pixels(matrix_b, sigma_b, threshold_min, threshold_max)
    matrix_c, sigma_c = process_bad_pixels(matrix_c, sigma_c, threshold_min, threshold_max)

    # Correct for beam stop, sum N files for non-normalized sensitivities
    matrix_d, sigma_matrix_d = calculate_pixel_wise_sensitivity(matrix_a, sigma_a,
                                                                matrix_b, sigma_b,
                                                                matrix_c, sigma_c)

    # Normalize pixel-wise sensitivities by weighting-average
    sensitivities, sensitivities_error, avg_sens, avg_sigma_sens = normalize_sensitivities(matrix_d,
                                                                                           sigma_matrix_d)

    # Compare with gold data
    verify_result(sensitivities, sensitivities_error)

    return


if __name__ == "__main__":
    pytest.main()
