import numpy as np
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
# from mantid.simpleapi import LoadEmptyInstrument, AddSampleLog
from drtsans.iq import bin_iq_into_linear_q1d, bin_iq_into_logarithm_q1d, IofQCalculator, BinningMethod
from drtsans.momentum_transfer_factory import calculate_q_dq

# This test implements issue #205 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
# DEV - Wenduo Zhou <petersonpf@ornl.gov> and Joe Osborn <osbornjd@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Lisa

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
                           [95, 97, 36 ]])

    l_matrix_c = np.array([[99, 400, 93],
                           [np.nan, 105, 94],
                           [95, 99, 36]])

    # Calculate uncertainties
    error_matrix_a = np.sqrt(l_matrix_a)
    error_matrix_b = np.sqrt(l_matrix_b)
    error_matrix_c = np.sqrt(l_matrix_c)

    return l_matrix_a, error_matrix_a, l_matrix_b, error_matrix_b, l_matrix_c, error_matrix_c


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
    # a = sum I(m, n) / sigma^2(m, n)
    weighted_sum = np.sum(normalized_data / normalized_error**2)
    # b = sum 1 / sigma^2(m, n)
    weights_square = np.sum(1. / normalized_error**2)
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

    return avg_norm_data, avg_norm_error,weighted_average, weighted_average_error


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


def calculate_scalar_sensitivity(d_matrix, sigma_d_matrix):
    """Do weighted average to pixel-wise sensitivities and propagate the error
    And then apply the average to sensitivity

    S_avg = \sum_{m, n}{D(m, n) / sigma^2(m, n)} / \sum_{m, n}{1 / sigma^2(m, n)}

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
    sens_avg = np.sum(d_matrix / sigma_d_matrix) / \
                      np.sum(1 / sigma_d_matrix)

    # Normalize pixel-wise sensitivities
    sensitivities = d_matrix / sens_avg

    # Calculate scalar sensitivity's error
    # sigma_S_avg = sqrt(1 / \sum_{m, n}(1 / sigma_D(m, n)^2))
    sigma_sens_avg = np.sqrt(1 / np.sum(1 / sigma_d_matrix))

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

    # Normalize by monitor
    matrix_a, sigma_a = test_data_set[0], test_data_set[1]
    matrix_a, sigma_a = normalize_by_monitor(matrix_a, sigma_a)

    # Find weighted average and error
    calculate_weighted_average_with_error()
    calculate_weighted_average_error()

    # Apply bad pixel threshold to the data
    process_bad_pixels(threshold_min, threshold_max)

    # Correct for beam stop
    single_file = calculate_pixel_wise_sensitivity()

    sensitivity = calculate_weighted_average_with_error()
    sen_error = calculate_weighted_average_error()

    # Compare with gold data
    verify_result(sensitivity, sen_error)

    return


if __name__ == "__main__":
    pytest.main()
