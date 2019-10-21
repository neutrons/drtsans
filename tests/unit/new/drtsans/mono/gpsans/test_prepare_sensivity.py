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


def normalize_by_monitor(data, error):
    """

    Parameters
    ----------
    data
    error

    Returns
    -------

    """

    return


def calculate_weighted_average():
    """Calculate weighted average

    weighted integration is same as iq.do_weighted_integration

    Returns
    -------
    ndarray
    """

    return


def calculate_weighted_average_error():
    """Propagate error for the weighted average

    Returns
    -------
    ndarray
    """

    return


def process_bad_pixels(threshold_min, threshold_max):
    """Process bad pixels

    If any pixel with counts falling out of allowed threshold, i.e., out of range (min, max)
    they will be specified as bad pixels.
    For any bad pixel, the counts will then be set to 'inf'

    Parameters
    ----------
    threshold_min: float
        minimum value of allowed value
    threshold_max: float
        maximum value of allowed value

    Returns
    -------

    """

    return

def correct_beam_stop():
    """Correct for beam stop by adding 3 files together

    Calculate Pixel-wise Average of 3 files to create the new summed file for
    doing the sensitivity correction

    Returns
    -------

    """
    return


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
    calculate_weighted_average()
    calculate_weighted_average_error()

    # Apply bad pixel threshold to the data
    process_bad_pixels(threshold_min, threshold_max)

    # Correct for beam stop
    single_file = correct_beam_stop()

    sensitivity = calculate_weighted_average()
    sen_error = calculate_weighted_average_error()

    # Compare with gold data
    verify_result(sensitivity, sen_error)

    return


if __name__ == "__main__":
    pytest.main()
