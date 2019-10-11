import numpy as np
import pytest


def assign_2d_bin_is(qx_array, qy_array, num_x_bins, num_y_bins):
    """Assign each pixel with a proper bin number

    Note: this method is exactly the same as Lisa's and is used to verify np.histogram2d()

    Note: 'is' in the method name stands for Instrument scientist (Lisa)

    Parameters
    ----------
    qx_array: ndarray (N, M)
        Qx
    qy_array: ndarray (N, M)
        Qy
    num_x_bins: integer
        number of bins on Qx array
    num_y_bins: integer
        number of bins on Qy array
    Returns
    -------

    """
    # Calculate delta(q_x) and delta(q_y)
    qx_bin_size = (np.max(qx_array) - np.min(qx_array)) / (num_x_bins - 1)
    qy_bin_size = (np.max(qy_array) - np.min(qy_array)) / (num_y_bins - 1)

    # Assign
    qx_index_array = np.ceil(qx_array / qx_bin_size + 0.5 * (1 + (np.max(qx_array) - np.min(qx_array)) / qx_bin_size))
    qy_index_array = np.ceil(qy_array / qy_bin_size + 0.5 * (1 + (np.max(qy_array) - np.min(qy_array)) / qy_bin_size))

    return qx_index_array, qy_index_array


def bin_iq_2d_is(i_array, qx_array, qy_array, num_x_bins, num_y_bins):
    """Bin I(Qx, Qy) for 2-dimensional Q

    Note: this method implementes Lisa's algorithm and is used to compare np.histogram2d()

    Parameters
    ----------
    i_array
    qx_array
    qy_array
    num_x_bins
    num_y_bins

    Returns
    -------

    """
    # Assign Qx and Qy
    qx_index_array, qy_index_array = assign_2d_bin_is(qx_array, qy_array, num_x_bins, num_y_bins)
    qx_index_array = qx_index_array.astype('int32')
    qy_index_array = qy_index_array.astype(int)

    # Bin (Qx', Qy')
    # TODO - continue from here
    raise NotImplementedError('ASAP as next')


if __name__ == "__main__":
    pytest.main()
