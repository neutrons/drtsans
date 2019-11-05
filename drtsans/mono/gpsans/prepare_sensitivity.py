"""
Module for algorithms to prepare sensitivity for instrument with moving detector
"""
import numpy as np

# Functions exposed to the general user (public) API
__all__ = ['prepare_sensitivity']


# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
def prepare_sensitivity(flood_data_matrix, flood_sigma_matrix, monitor_counts, threshold_min, threshold_max):
    """Prepare sensitivity for moving detector

    Data files are processed such that intensities and errors are stored in numpy.ndarray with shape (N, M), where
    - N: number of data files to calculate sensitivities
    - M: number of pixels (aka spectra) in instrument's detector;
         The 2D data from 2D detector are flattened to 1D in implementation

    Prerequisite of the input files:
    - top and bottom of the detector shall be masked (set to value as NaN) due to edge effects
    - in each file, beam center shall be found  and masked out
    - errors are then calculated from the flood intensities

    Workflow to prepare sensitivities
    - normalize the flood field data by monitor
    - find weighted average for each fie and error
    - apply bad pixel threshold to each file
    - correct for beam stop and add all the flood files together to non-normalized sensitivities
    - apply weighted average to sensitivities

    Parameters
    ----------
    flood_data_matrix : ~numpy.ndaray
        multiple set of flood data intensities with shape = N, M
    flood_sigma_matrix : ~numpy.ndaray
        multiple set of flood data intensities' error with shape = N, M
    monitor_counts : ~numpy.ndaray
        monitor counts for each data file in 1D array with size N
    threshold_min : float
        minimum allowed detector counts to mask out 'bad' pixels
    threshold_max : float
        maximum allowed detector counts to mask out 'bad' pixels

    Returns
    -------
    ~numpy.ndaray, ~numpy.ndaray
        sensitivities, sensitivities error

    """
    # normalize the flood field data by monitor
    # inputs: (N, M) array; outputs: (N, M) array
    flood_data_matrix, flood_sigma_matrix = _normalize_by_monitor(flood_data_matrix, flood_sigma_matrix,
                                                                  monitor_counts)

    # find weighted average for each fie and error
    # inputs: (N, M) array; outputs: (N, M) array
    returns = _calculate_weighted_average_with_error(flood_data_matrix, flood_sigma_matrix)
    flood_data_matrix = returns[0]
    flood_sigma_matrix = returns[1]

    # apply bad pixel threshold to each file
    # inputs: (N, M) array; outputs: (N, M) array
    flood_data_matrix, flood_sigma_matrix = _process_bad_pixels(flood_data_matrix, flood_sigma_matrix,
                                                                threshold_min, threshold_max)

    # correct for beam stop and add all the flood files together to non-normalized sensitivities
    raw_sensitivities, raw_sensitivities_error = _calculate_pixel_wise_sensitivity(flood_data_matrix,
                                                                                   flood_sigma_matrix)
    #
    # # apply weighted average to sensitivities
    # sensitivities, sensitivities_error,  sens_avg, sigma_sens_avg = _normalize_sensitivities(raw_sensitivities,
    #                                                                                          raw_sensitivities_error)
    # print('[DEBUG] Sensitivity Avg = {}, Sigma Avg = {}'.format(sens_avg, sigma_sens_avg))
    # return sensitivities, sensitivities_error

    return raw_sensitivities, raw_sensitivities_error


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
    assert monitor_counts.shape == (flood_data.shape[0], 1), 'Monitor counts must be in shape as ({}, 1} ' \
                                                             'but not {}'.format(flood_data.shape[0],
                                                                                 monitor_counts.shape)

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
    data_no_nan = normalized_data[~np.isnan(normalized_data)]
    print('[DEBUG] Data with NaN removed: {}\n....... shape = {}'.format(data_no_nan, data_no_nan.shape))
    weighted_sum = np.nansum(normalized_data / normalized_error**2, axis=1)  # summing in row
    # b = sum 1 / sigma^2(m, n)
    weights_square = np.nansum(1. / normalized_error**2, axis=1)
    # Avg = a / b
    weighted_average = weighted_sum / weights_square
    weighted_average = weighted_average.reshape((normalized_data.shape[0], 1))  # reshape to (N, 1) for division
    # sigma Avg = 1 / sqrt(b)
    weighted_average_error = 1. / np.sqrt(weights_square)
    weighted_average_error = weighted_average_error.reshape((normalized_data.shape[0], 1))
    print('[DEBUG] Average = {}, Error = {}'.format(weighted_average, weighted_average_error))

    # Normalize data by weighted-average
    avg_norm_data = normalized_data / weighted_average

    # Propagate uncertainties: sigma S(m, n) = I(m, n) / avg * [(error(m, n)/I(m, n))^2 + (sigma Avg/Avg)^2]^1/2
    # in the sqrt operation, first term is a N x M array and second term is a N x 1 array
    avg_norm_error = normalized_data / weighted_average * np.sqrt((normalized_error / normalized_data)**2
                                                                  + (weighted_average_error / weighted_average)**2)

    assert normalized_data.shape == avg_norm_error.shape

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


def _calculate_pixel_wise_sensitivity(flood_data, flood_error):
    """Calculate pixel-wise average of N files to create the new summed file for doing sensitivity correction

    # data_a, data_a_error, data_b, data_b_error, data_c, data_c_error
    D(m, n) = A_F(m, n) + B_F(m, n) + C_F(m, n) with average weight

    Calculate Pixel-wise Average of 3 files to create the new summed file for
    doing the sensitivity correction

    Parameters
    ----------
    flood_data : ~numpy.ndarray
        processed multiple flood files in an N x M array
    flood_error : ~numpy.ndarray
        processed multiple flood files' error in an N x M array

    Returns
    -------
    ~numpy.nparray, ~numpy.nparray
        non-normalized sensitivities, non-normalized sensitivities error
        1D array as all the flood files are summed

    """
    # Create sensitivities and sigmas matrices
    sensitivities = np.zeros_like(flood_data[0])
    sensitivities_error = np.zeros_like(flood_data[0])

    # Calculate D'(i, j)    = sum_{k}^{A, B, C}M_k(i, j)/s_k^2(i, j)
    #           1/s^2(i, j) = sum_{k}^{A, B, C}1/s_k^2(i, j)
    # flood_data.shape[0]: number of flood files
    # flood_data.shape[1]: number of pixels
    for ipixel in range(flood_data.shape[1]):
        # For each pixel: sum up along axis = 1
        d_ij_arr = flood_data[:, ipixel]
        s_ij_arr = flood_error[:, ipixel]

        print('[DEBUG] Pixel {}: D_ij = {}'.format(ipixel, d_ij_arr))
        print('[DEBUG] INF: {}'.format(np.where(np.isinf(d_ij_arr))))

        if len(np.where(np.isinf(d_ij_arr))) > 0:
            # In case there is at least an inf in this subset of data, set sensitivities to -inf
            sensitivities[ipixel] = -np.inf
            sensitivities_error[ipixel] = -np.inf
        else:
            # Do weighted summation to the subset by excluding the NaN
            s_ij = np.nansum(1. / s_ij_arr ** 2)
            d_ij = np.nansum(d_ij_arr / s_ij_arr ** 2) / s_ij
            s_ij = np.sqrt(s_ij)
            sensitivities[ipixel] = d_ij
            sensitivities_error[ipixel] = 1. / s_ij
        # END-IF-ELSE
    # END-FOR

    return sensitivities, sensitivities_error


def _normalize_sensitivities(d_matrix, sigma_d_matrix):
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
