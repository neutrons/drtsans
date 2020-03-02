import numpy as np
from collections import namedtuple


# Define structure (namedtuple) for binning parameters: min, max, number of bins
# bins shall be integer as number of bins
BinningParams = namedtuple('BinningParams', 'min max bins')
# Define structure (namedtuple) for bins: bin edges and bin boundaries
# Both bin edge and bin boundaries shall be 1-dimensional 1D array and 'edges' size is 1 larger than centers
Bins = namedtuple('Bins', 'edges centers')


def determine_1d_linear_bins(x_min, x_max, bins):
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
    # Check input x min and x max
    if x_min is None or x_max is None or x_min >= x_max:
        raise RuntimeError('x min {} and x max {} must not be None and x min shall be less than x max'
                           ''.format(x_min, x_max))

    # Calculate Q step size
    delta_x = (x_max - x_min) / bins
    # Determine bin edges
    bin_edges = np.arange(bins + 1).astype('float') * delta_x + x_min
    # Determine bin centers from edges
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

    # Construct Bins instance
    linear_bins = Bins(bin_edges, bin_centers)

    return linear_bins


def determine_1d_log_bins(x_min, x_max, n_bins_per_decade=None, n_bins=None,
                          decade_on_center=False, even_decade=False):
    """Determine logarithm bins

    Including bin edge and bin center.

    n_bins_per_decade and n_bins cannot be defined simultaneously, but one of them must be specified.

    Parameters
    ----------
    x_min : float
        Positive float for minimum
    x_max : float
        Positive float
    n_bins_per_decade : int or None
        Positive integer for number of step per decade. Total number of bins will be this value multiplied by
        number of decades from X min to X max
    n_bins : int or None
        Positive integer for total number of bins.
    decade_on_center : bool
        Flag to have the min X and max X on bin center; Otherwise, they will be on bin boundary
    even_decade : bool
        Flag to have even decade for minimum and maximum value in the generated bins

    Returns
    -------
    ~drtsans.iq.Bins
        Bins including bin centers and bin edges

    """
    # Calculate C min and max on decade and contain X min and X max in the range
    # From Chapter 11 of the master document
    if even_decade:
        # both c_min and c_max shall be on even decade
        # x_min and x_max shall be enclosed in the range of (10^c_min, 10^c_max)
        c_min = np.log10(x_min)
        c_max = np.ceil(np.log10(x_max))
    else:
        # c_min and c_max are from x_min and x_max directly
        c_min = np.log10(x_min)
        c_max = np.log10(x_max)

    # Calculate total number of bins
    if n_bins_per_decade is not None and n_bins_per_decade > 0 and n_bins is None:
        # Number of bins shall contain all the decade from c_min to c_max
        # Thus, floor to c_min and ciel to c_max shall make sure the calculation is correct
        total_num_bins = (int(np.ceil(c_max) - np.floor(c_min))) * n_bins_per_decade
    elif n_bins_per_decade is None and n_bins is not None and n_bins > 0:
        # Use user specified total number of bins
        total_num_bins = n_bins
    else:
        # Non-supported case
        raise RuntimeError('n_bins_per_decade ({}) and n_bins ({}) cannot be both specified or not specified.'
                           'and they must be positive integers'.format(n_bins_per_decade, n_bins))

    # Calculate Delta L: Equation 11.28 (master document)
    delta_l = (c_max - c_min) / total_num_bins

    # Define an array of k, i.e., [0, 1, 2, ...]
    if decade_on_center:
        # x_min and 10^{c_max} on the first and last bin center
        # Equation 11.31: number of bins will be N + 1 for bin on the center
        total_num_bins += 1
        vec_k = np.arange(total_num_bins).astype(float)
        bin_centers = np.power(10, delta_l * vec_k + c_min)
    else:
        # x_min and 10^{c_max} on the first and last bin boundary
        # Equation 11.29
        vec_k = np.arange(total_num_bins).astype(float)
        bin_centers = np.power(10, delta_l * (vec_k + 0.5) + c_min)
    # END-IF-ELSE

    # Calculate bin boundaries from bin center
    # Equation 11.30
    bin_edges = np.ndarray(shape=(bin_centers.shape[0] + 1,), dtype=float)
    bin_edges[1:-1] = 0.5 * (bin_centers[:-1] + bin_centers[1:])

    # Set the min and max of bins explicitly
    if decade_on_center:
        # x_min and 10^{c_max} are on the first and last bin center
        # then first and last bin edges/boundaries are defined as
        # 10^{C_min - deltaL / 2 } and 10^{C_max + deltaL / 2}
        # according to the paragraph after equation 11.31
        bin_edges[0] = np.power(10, c_min - 0.5 * delta_l)
        bin_edges[-1] = np.power(10, c_max + 0.5 * delta_l)
    elif even_decade:
        # x_min and 10^{c_max} on the first and last bin boundary
        # then first and last bin edges/boundaries are defined as
        # Q_min and Q_max (or X_min and X_max in generalized form)
        bin_edges[0] = 10**c_min
        bin_edges[-1] = 10**c_max
    else:
        # use user x min and x max
        bin_edges[0] = x_min
        bin_edges[-1] = x_max
    # END-IF-ELSE

    # Construct Bins instance
    log_bins = Bins(bin_edges, bin_centers)

    return log_bins
