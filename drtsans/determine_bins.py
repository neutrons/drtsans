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

    Returns
    -------
    ~drtsans.iq.Bins
        Bins including bin centers and bin edges

    """
    # Check input x min and x max
    if x_min is None or x_max is None or x_min >= x_max:
        raise RuntimeError('x min {} and x max {} must not be None and x min shall be less than x max'
                           ''.format(x_min, x_max))
    # force the number of bins to be an integer and error check it
    bins = int(bins)
    if bins <= 0:
        raise ValueError('Encountered illegal number of bins: {}'.format(bins))

    # Calculate Q step size
    delta_x = float((x_max - x_min) / bins)
    # Determine bin edges
    bin_edges = np.arange(bins + 1).astype('float') * delta_x + x_min
    # Determine bin centers from edges
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

    # Construct Bins instance
    linear_bins = Bins(bin_edges, bin_centers)

    return linear_bins


def determine_1d_log_bins_new(x_min, x_max, decade_on_center,
                              n_bins_per_decade=None, n_bins=None):
    """

    Parameters
    ----------
    x_min: float
        minimum value of X in the bins
    x_max: float
        maximum value of X in the bins
    decade_on_center: bool
        flag that data must be centered on decades
    n_bins_per_decade: int, None
        density of points (number of data points per decade)
    n_bins: int, None
        total number of points in the output

    Returns
    -------

    """
    # Check inputs
    if n_bins_per_decade is None and n_bins is None:
        raise RuntimeError('blabla')
    elif n_bins_per_decade is not None and n_bins is not None:
        raise RuntimeError('blabla')
    # only allow either n_bins or n_bins_per_decade

    # Calculate Q min, number of total bins and number of steps
    if n_bins_per_decade is not None:
        # user specifies number of bins per decade

        # determine X min
        if decade_on_center:
            x_ref = _calculate_x_ref(x_min, n_bins_per_decade)
            print(f'reference X = {x_ref}')
        else:
            x_ref = x_min

        # calculate step size
        n_step = 10**(1 / n_bins_per_decade)
        print(f'n_step = {n_step}')

        # calculate number of bins
        n_bins = _calculate_n_bins(x_ref, x_max, n_step, n_bins_per_decade)
        print(f'total number of bins = {n_bins}')

    else:
        # user specifies number of total bins

        # case that is not supported
        if decade_on_center:
            assert n_bins_per_decade is not None, 'For option decade_on_center, number of bins per decade ' \
                                                  'is required'
        x_ref = x_min

        # calculate bin step size
        # Equation 11.33
        n_step = 10**((np.log10(x_max / x_ref)) / (n_bins - 1))
        print(f'n_steps = {n_step}')

    # Calculate kay
    kay = (n_step - 1) / (n_step + 1)

    # Calculate bin centers
    # init an array ranging from 0 to (n_bins - 1)
    bin_centers = np.arange(n_bins).astype('float64')
    # Equation 11.34: Q_k = Q_ref * 10^(k * delta L)
    bin_centers = x_ref * n_step**bin_centers

    # Calculate bin edges (aka boundaries)
    # Equation 11.35
    bin_edges = np.ndarray(shape=(n_bins + 1, ), dtype='float64')
    # calculate left edges (i.e., right edges except last one), i.e., Q_{k-1, max} = Q_{k, min}
    bin_edges[:-1] = bin_centers[:] - kay * bin_centers[:]
    # calculate last right edge
    bin_edges[-1] = bin_centers[-1] + kay * bin_centers[-1]

    # Form output as Bins object
    log_bins = Bins(bin_edges, bin_centers)

    return log_bins


def _calculate_x_ref(x_min, n_bins_per_decade):
    """Calculate reference X (minimum X) by Equation in Chapter 11
    x_ref = 10^((1/N_bins_decade)*(round(N_bins_decade*log10(Q_min))))
    Parameters
    ----------
    x_min: float
        minimum x value
    n_bins_per_decade: int
        point density

    Returns
    -------

    """
    ref_x = 10 ** ((1. / n_bins_per_decade) * (np.round(n_bins_per_decade * np.log10(x_min))))

    return ref_x


def _calculate_n_bins(x_min, x_max, n_step, bin_density):
    """Calculate number of total bins by implementing
    Equation 11.32

    N_bins = floor(ceil((log10(Q_max/Calc_Q_min) + log10((N_step+1)/2.0))/log10(N_step)))

    Parameters
    ----------
    x_min
    x_max
    n_step

    Returns
    -------

    """
    n_bins = np.floor(np.ceil((np.log10(x_max / x_min) + np.log10((n_step + 1) * 0.5)) / np.log10(n_step)))

    """
        t0 = np.log10(x_max / x_min)
        t1 = bin_density * np.log10((1 + 10**(1. / bin_density)) / 2.0)
        print(f't0 = {t0}')
        print(f't1 = {t1}')
        n_bins_x = np.ceil(t0 + t1)
        print(f'[DEBUG] total bins = {n_bins}  vs {n_bins_x}')
    """

    # to avoid round off error such that n_bins = |n_bins| + epsilon, where epsilon is an infinitesimally
    # small value
    n_bins = int(n_bins + 1E-5)

    return n_bins


def determine_1d_log_bins_deprecated(x_min, x_max, n_bins_per_decade=None, n_bins=None,
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
        c_min = np.floor(np.log10(x_min))
        c_max = np.ceil(np.log10(x_max))
    else:
        # c_min and c_max are from x_min and x_max directly
        c_min = np.log10(x_min)
        c_max = np.log10(x_max)

    # Calculate total number of bins
    if n_bins_per_decade is not None and n_bins_per_decade > 0 and n_bins is None:
        # Number of bins shall contain all the decade from c_min to c_max
        # Thus, floor to c_min and ciel to c_max shall make sure the calculation is correct
        total_num_bins = int(np.ceil((c_max - c_min) * n_bins_per_decade))
        delta_l = 1./n_bins_per_decade
    elif n_bins_per_decade is None and n_bins is not None and n_bins > 0:
        # Use user specified total number of bins
        total_num_bins = n_bins
        # Calculate Delta L: Equation 11.28 (master document)
        delta_l = (c_max - c_min) / total_num_bins
    else:
        # Non-supported case
        raise RuntimeError('n_bins_per_decade ({}) and n_bins ({}) cannot be both specified or not specified.'
                           'and they must be positive integers'.format(n_bins_per_decade, n_bins))

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

    # Calculate bin boundaries (edges) from bin center
    # Equation 11.30 revised by Ken
    n_step = 10 ** (1 / n_bins_per_decade)
    kay = (n_step - 1) / (n_step + 1)

    bin_edges = np.ndarray(shape=(bin_centers.shape[0] + 1,), dtype=float)
    bin_edges[:-1] = bin_centers[:] - kay * bin_centers[:]

    # Set the min and max of bins explicitly
    if decade_on_center:
        # x_min and 10^{c_max} are on the first and last bin center
        # then first and last bin edges/boundaries are defined as
        # 10^{C_min - deltaL / 2 } and 10^{C_max + deltaL / 2}
        # according to the paragraph after equation 11.31
        # bin_edges[0] = np.power(10, c_min - 0.5 * delta_l)
        bin_edges[-1] = bin_centers[-1] + kay * bin_centers[-1]
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

    # recalculate the last center
    if bin_centers[-1] >= bin_edges[-1]:
        bin_centers[-1] = 0.5 * (bin_edges[-1] + bin_edges[-2])

    # Construct Bins instance
    log_bins = Bins(bin_edges, bin_centers)

    return log_bins
