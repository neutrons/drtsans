import numpy as np
from collections import namedtuple


# Define structure (namedtuple) for binning parameters: min, max, number of bins
# bins shall be integer as number of bins
BinningParams = namedtuple('BinningParams', 'min max bins')
# Define structure (namedtuple) for bins: bin edges and bin boundaries
# Both bin edge and bin boundaries shall be 1-dimensional 1D array and 'edges' size is 1 larger than centers
Bins = namedtuple('Bins', 'edges centers')


def determine_1d_linear_bins(x_min, x_max, bins, x_array=None):
    """Determine linear bin edges and centers

    Parameters
    ----------
    x_min : float
        Q min of bin edge
    x_max : float or None
        Q max of bin edge
    bins : integer
        number of bins
    x_array: numpy.ndarray
        1-dimensional numpy array for X-axis such as Q, Qx, Qy or annular angles
        The default is None if min X and max X are specified by users

    Returns
    -------
    ~drtsans.iq.Bins
        Bins including bin centers and bin edges

    """
    # Define q min and q max if they are not specified
    try:
        x_min = x_array.min() if x_min is None else x_min
        x_max = x_array.max() if x_max is None else x_max
    except AttributeError as att_error:
        raise RuntimeError('X Min ({}) or X Max ({}) must be specified if X-array is not given: FYI {}'
                           ''.format(x_min, x_max, att_error))

    # Calculate Q step size
    delta_x = (x_max - x_min) / bins
    # Determine bin edges
    bin_edges = np.arange(bins + 1).astype('float') * delta_x + x_min
    # Determine bin centers from edges
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

    # Construct Bins instance
    linear_bins = Bins(bin_edges, bin_centers)

    return linear_bins


def determine_1d_log_bins(x_min, x_max, step_per_decade):
    """Determine the logarithm bins

    The algorithm to calculate bins are from 11.28 and 11.29 in master document

    Parameters
    ----------
    x_min
    x_max
    step_per_decade: float
        step per decade (ex. 0.1 to 1.0 is one decade); denoted as 'j' in document
    Returns
    -------
    Bins
        bin centers, bin edges

    """
    # C_max = 10^ceil(log{Q_max}): nearest equal or larger power of 10 for C_max
    c_max = 10 ** (np.ceil(np.log10(x_max)))
    # C_min = 10^floor(log(Q_min)): nearest equal or smaller power of 10 for C_min
    c_min = 10 ** (np.floor(np.log10(x_min)))
    # log step
    delta_l = 1. / step_per_decade
    # number of total bins
    num_bins = int(np.log10(c_max / c_min)) * step_per_decade

    # Determine Q bin centers
    bin_centers = np.arange(num_bins)
    bin_centers = bin_centers.astype(float)
    # bin centers: 10^{delta_l * (k + 0.5)} * c_min
    bin_centers = 10 ** (delta_l * (bin_centers + 0.5)) * c_min
    # Determine Q bin edges
    bin_edges = np.zeros((bin_centers.shape[0] + 1), float)
    bin_edges[0] = c_min
    bin_edges[1:-1] = 0.5 * (bin_centers[:-1] + bin_centers[1:])
    bin_edges[-1] = c_max

    # Initialize Bins instance
    log_bins = Bins(bin_edges, bin_centers)

    return log_bins

