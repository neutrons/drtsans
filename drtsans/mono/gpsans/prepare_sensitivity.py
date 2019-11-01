"""
Module for algorithms to prepare sensitivity for instrument with moving detector
"""
import numpy as np

# Functions exposed to the general user (public) API
__all__ = ['prepare_sensitivity']


def prepare_sensitivity(flood_data_matrix, flood_sigma_matrix, monitor_counts, threshold_min, threshold_max):
    """Prepare sensitivity for moving detector

    Parameters
    ----------
    flood_data : ndarray
        multiple set of flood data.  flood_data.shape = (M, N) as M set of N pixels

    Returns
    -------

    """
    # Normalize by
    flood_data_matrix, flood_sigma_matrix = _normalize_by_monitor(flood_data_matrix, flood_sigma_matrix,
                                                                  monitor_counts)

    _calculate_weighted_average_with_error(flood_data_matrix, flood_sigma_matrix)

    _process_bad_pixels(flood_data_matrix, flood_sigma_matrix, threshold_min, threshold_max)

    return


def _normalize_by_monitor(flood_data, flood_data_error, monitor_counts):
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
    # Check monitor counts shape and convert if necessary
    if len(monitor_counts.shape) == 1:
        monitor_counts = monitor_counts.reshape((len(monitor_counts), 1))
    assert monitor_counts.shape == flood_data.shape[0], 1

    return flood_data / monitor_counts, flood_data_error / monitor_counts


def _calculate_weighted_average_with_error(normalized_data, normalized_error):
    """Calculated weighted average for normalized flood data and error

    Average = (sum I(m, n) / sigma^(m, n)) / sum 1 / sigma^2

    Parameters
    ----------
    normalized_data : ndarray
        normalized flood data
        multiple set of flood data.  flood_data.shape = (M, N) as M set of N pixels
    normalized_error : ndarray
        normalized flood data's error
        multiple set of flood data.  flood_data.shape = (M, N) as M set of N pixels

    Returns
    -------
    ndarray, ndarray, float, float
        data normalized by average, data's error normalized by average, Average, sigma(Average)

    """
    # Calculate weighted average
    # a = sum_{m, n} I(m, n) / sigma^2(m, n)
    weighted_sum = np.sum(normalized_data[~np.isnan(normalized_data)] /
                          normalized_error[~np.isnan(normalized_error)]**2,
                          axis=1)
    # b = sum 1 / sigma^2(m, n)
    weights_square = np.sum(1. / normalized_error[~np.isnan(normalized_error)]**2,
                            axis=1)
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


def _process_bad_pixels(data, data_error, threshold_min, threshold_max):
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
