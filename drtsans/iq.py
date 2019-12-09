from collections import namedtuple
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/dataobjects.py
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/docs/drtsans/dataobjects.rst
from drtsans.dataobjects import IQazimuthal, IQmod
from enum import Enum
import numpy as np
# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')

__all__ = ['bin_intensity_into_q1d', 'select_i_of_q_by_wedge',
           'bin_annular_into_q1d', 'bin_intensity_into_q2d', 'BinningMethod', 'check_iq_for_binning',
           'determine_1d_linear_bins', 'determine_1d_log_bins', 'BinningParams']


# Define structure (namedtuple) for binning parameters: min, max, number of bins
# bins shall be integer as number of bins
BinningParams = namedtuple('BinningParams', 'min max bins')
# Define structure (namedtuple) for bins: bin edges and bin boundaries
# Both bin edge and bin boundaries shall be 1-dimensional 1D array and 'edges' size is 1 larger than centers
Bins = namedtuple('Bins', 'edges centers')


class BinningMethod(Enum):
    """
    Binning method
    """
    NOWEIGHT = 1   # no-weight binning
    WEIGHTED = 2   # weighted binning


def check_iq_for_binning(i_of_q):
    """Check I(Q) for binning.

    Binning I(Q) assumes that
    1. there is no NaN or Infinity in intensities
    2. there is no NaN, Infinity or Zero in intensity errors

    :exception : RuntimeError
        raise exception if input I(Q) does not meet assumption

    :param i_of_q:  ~drtsans.dataobjects.IQmod or IQazimuthal
        I(Q)

    """
    error_message = ''

    # Check intensity
    if np.where(np.isnan(i_of_q.intensity))[0].size > 0:
        error_message += 'Intensity has NaN: {}\n'.format(np.where(np.isnan(i_of_q.intensity))[0])
    if np.where(np.isinf(i_of_q.intensity))[0].size > 0:
        error_message += 'Intensity has Inf: {}\n'.format(np.where(np.isnan(i_of_q.intensity))[0])

    # Check error
    if np.where(np.isnan(i_of_q.error))[0].size > 0:
        error_message += 'Intensity error has NaN: {}\n'.format(np.where(np.isnan(i_of_q.error))[0])
    if np.where(np.isinf(i_of_q.error))[0].size > 0:
        error_message += 'Intensity error has Inf: {}\n'.format(np.where(np.isnan(i_of_q.error))[0])
    if np.where(np.abs(i_of_q.error) < 1E-20)[0].size > 0:
        error_message += 'Intensity error has zero {}\n'.format(np.where(np.abs(i_of_q.error) < 1E-20)[0])

    if len(error_message) > 0:
        raise RuntimeError('Input I(Q) for binning does not meet assumption:\n{}'.format(error_message))


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


def bin_intensity_into_q1d(i_of_q, q_bins, bin_method=BinningMethod.WEIGHTED):
    """Binning I(Q) from scalar Q (1D) with linear binning on Q

    Replace intensity, intensity_error, scalar_q, scalar_dq by IQmod
    Replace bins, q_min=None, q_max=None by BinningParams
    bins: number of bins for linear binning; step per decade for logarithm binning
    q_min : Default to min(scalar_q)
    q_max : Default to max(scalar_q)

    Parameters
    ----------
    i_of_q : ~drtsans.dataobjects.IQmod
        Scalar I(Q) including intensity, intensity_error, scalar_q, scalar_dq in 1d nparray
        including: intensity error mod_q delta_mod_q
    q_bins : Bins
        namedtuple for arbitrary bin edges and bin centers
    bin_method : ~drtsans.BinningMethod
        weighted binning or no-weight binning method

    Returns
    -------
    drtsans.dataobjects.IQmod
        the one dimensional data as a named tuple
    """
    # Check input I(Q) whether it meets assumptions
    check_iq_for_binning(i_of_q)

    # bin I(Q)
    if bin_method == BinningMethod.WEIGHTED:
        # weighed binning
        binned_q = _do_1d_weighted_binning(i_of_q.mod_q, i_of_q.delta_mod_q, i_of_q.intensity, i_of_q.error,
                                           q_bins.centers, q_bins.edges)
    else:
        # no-weight binning
        binned_q = _do_1d_no_weight_binning(i_of_q.mod_q, i_of_q.delta_mod_q, i_of_q.intensity, i_of_q.error,
                                            q_bins.centers, q_bins.edges)

    return binned_q


def select_i_of_q_by_wedge(i_of_q, min_wedge_angle, max_wedge_angle):
    """Select a sub set of I(Q) by 2D wedge

    Parameters
    ----------
    i_of_q : ~drtsans.dataobjects.IQazimuthal
         "intensity": intensity, "error": sigma(I), "qx": qx, "qy": qy, "delta_qx": dqx, "delta_qy", dqy
    min_wedge_angle : float
        minimum value of theta/azimuthal angle for wedge
    max_wedge_angle : float
        maximum value of theta/azimuthal angle for wedge

    Returns
    -------
    drtsans.dataobjects.IQazimuthal
        subset of input I(Qx, Qy) with (Qx, Qy) inside defined wedge

    """
    # Calculate wedge angles for each I(Qx, Qy)
    # calculate azimuthal angles from -180 to 180 degrees
    azimuthal_array = np.arctan2(i_of_q.qy, i_of_q.qx) * 180. / np.pi
    # correct azimuthal angles to -90 to 270 degrees
    azimuthal_array[azimuthal_array < -90.] += 360.

    # Define the filter (mask/ROI) for pixels falling into preferred wedge
    wedge_indexes = (azimuthal_array > min_wedge_angle) & (azimuthal_array < max_wedge_angle)

    # Create a new IQazimuthal with selected subset
    wedge_i_of_q = IQazimuthal(i_of_q.intensity[wedge_indexes],
                               i_of_q.error[wedge_indexes],
                               i_of_q.qx[wedge_indexes],
                               i_of_q.qy[wedge_indexes],
                               i_of_q.delta_qx[wedge_indexes],
                               i_of_q.delta_qy[wedge_indexes])

    return wedge_i_of_q


def bin_annular_into_q1d(i_of_q, theta_bin_params, q_min=0.001, q_max=0.4, method=BinningMethod.NOWEIGHT):
    """Annular 1D binning

    Calculates: I(Q), sigma I and dQ by assigning pixels to proper azimuthal angle bins
    Given I(Qx, Qy) and will convert to I(Q) in the code

    Parameters
    ----------
    i_of_q :  drtsans.dataobjects.IQazimuthal
        I(Qx, Qy), sigma I(Qx, Qy), Qx, Qy, dQx and dQy
    theta_bin_params : ~drtsans.BinningParams
        binning parameters on annular angle 'theta'

        theta_min : float
            minimum value of theta/azimuthal angle
        theta_max : float
            maximum value of theta/azimuthal angle
        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in
            the given range (10 by default).  If `bins` is a sequence, it
            defines the bin edges, including the rightmost edge, allowing for
            non-uniform bin widths.  Values in `x` that are smaller than lowest
            bin edge areassigned to bin number 0, values beyond the highest bin
            are assigned to ``bins[-1]``.  If the bin edges are specified,
            the number of bins will be, (nx = len(bins)-1).

    q_min : float, optional
        , by default
    q_max : float, optional
        , by default
    method : ~drtsans.BinningMethod
        binning method, no-weight or weighed

    Returns
    -------
    drtsans.dataobjects.IQmod
        Annular-binned I(Q) in 1D

    """
    # Determine azimuthal angle bins (i.e., theta bins)
    theta_bins = determine_1d_linear_bins(theta_bin_params.min, theta_bin_params.max, theta_bin_params.bins)

    # Calculate theta array
    theta_array = np.arctan2(i_of_q.qy, i_of_q.qx) * 180. / np.pi
    # convert -0 to -180 to 180 to 360
    theta_array[np.where(theta_array < 0)] += 360.

    # Calculate Q from Qx and Qy
    q_array = np.sqrt(i_of_q.qx ** 2 + i_of_q.qy ** 2)
    # calculate dQ from dQx and dQy
    dq_array = np.sqrt(i_of_q.delta_qx ** 2 + i_of_q.delta_qy ** 2)

    # Filter by q_min and q_max
    allowed_q_index = (q_array > q_min) & (q_array < q_max)

    # binning
    if method == BinningMethod.NOWEIGHT:
        # no weight binning
        binned_iq = _do_1d_no_weight_binning(theta_array[allowed_q_index],
                                             dq_array[allowed_q_index],
                                             i_of_q.intensity[allowed_q_index],
                                             i_of_q.error[allowed_q_index],
                                             theta_bins.centers, theta_bins.edges)
    elif method == BinningMethod.WEIGHTED:
        # weighted binning
        # Check input I(Q) whether it meets assumptions
        check_iq_for_binning(i_of_q)
        binned_iq = _do_1d_weighted_binning(theta_array[allowed_q_index],
                                            dq_array[allowed_q_index],
                                            i_of_q.intensity[allowed_q_index],
                                            i_of_q.error[allowed_q_index],
                                            theta_bins.centers, theta_bins.edges)
    else:
        # not supported case
        raise RuntimeError('Binning method {} is not recognized'.format(method))

    return binned_iq


def _determine_1d_log_bins_lisa(q_min, q_max, step_per_decade):
    """Lisa's algorithm tot determine 1D log bins

    Parameters
    ----------
    q_min
    q_max
    step_per_decade: float
        step per decade (ex. 0.1 to 1.0 is one decade); denoted as 'j' in document
    Returns
    -------
    ndarray, ndarray
        bin centers, bin edges

    """
    # Calculate step and align q_min to q0, a decade (power of 10) nearest to q_min but less than q_min
    # 20191016 IS: "0.2% error is allowed.  This formula ensure that the number of steps per decade is respected"
    delta = np.power(10., 1. / step_per_decade)
    q0 = np.power(delta, np.floor(step_per_decade * np.log10(q_min)))

    # Determine number of bins
    num_bins = 1 + int(np.ceil(step_per_decade * np.log(q_max / q0) / np.log(10)))

    # Calculate bin centers
    bin_centers = np.arange(num_bins).astype('float')
    bin_centers = q0 * np.power(delta, bin_centers)

    # Calculate bin boundaries
    delta_q_array = 2. * (delta - 1) / (delta + 1) * bin_centers
    bin_edges = np.zeros((num_bins + 1,), dtype='float')
    bin_edges[1:] = bin_centers[:] + 0.5 * delta_q_array[:]
    bin_edges[0] = bin_centers[0] - 0.5 * delta_q_array[0]

    return bin_centers, bin_edges


def _do_1d_no_weight_binning(q_array, dq_array, iq_array, sigmaq_array, bin_centers, bin_edges):
    """ Bin I(Q) by given bin edges and do no-weight binning
    This method implements equation 11.34, 11.35 and 11.36 in master document.
    Parameters
    ----------
    q_array: ndarray
        scalar momentum transfer Q in 1D array flattened from 2D detector
    dq_array: ndarray
        scalar momentum transfer (Q) resolution in 1D array flattened from 2D detector
    iq_array: ndarray
        I(Q) in 1D array flattened from 2D detector
    sigmaq_array: ndarray
        sigma I(Q) in 1D array flattened from 2D detector
    bin_centers: numpy.ndarray
        bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
    bin_edges: numpy.ndarray
        bin edges
    Returns
    -------
    ~drtsans.dataobjects.IQmod
        IQmod is a class for holding 1D binned data.

    """
    # check input
    assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

    # Count number of Q in 'q_array' in each Q-bin when they are binned (histogram) to 'bin_edges'
    num_pt_array, bin_x = np.histogram(q_array, bins=bin_edges)

    # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
    i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array)

    # Square of summed uncertainties for each bin
    sigma_sqr_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=sigmaq_array ** 2)

    # Final I(Q):     I_k       = \frac{I_{k, raw}}{N_k}
    i_final_array = i_raw_array / num_pt_array
    # Final sigma(Q): sigmaI_k  = \frac{sigmaI_{k, raw}}{N_k}
    sigma_final_array = np.sqrt(sigma_sqr_array) / num_pt_array

    # Calculate Q resolution of binned
    binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
    bin_q_resolution = binned_dq / num_pt_array

    # Get the final result by constructing an IQmod object defined in ~drtsans.dataobjects.
    # IQmod is a class for holding 1D binned data.
    return IQmod(intensity=i_final_array, error=sigma_final_array,
                 mod_q=bin_centers, delta_mod_q=bin_q_resolution)


def _do_1d_weighted_binning(q_array, dq_array, iq_array, sigma_iq_array, bin_centers, bin_edges):
    """ Bin I(Q) by given bin edges and do weighted binning

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

    Parameters
    ----------
    q_array: ndarray
        scalar momentum transfer Q in 1D array flattened from 2D detector
    dq_array: ndarray
        scalar momentum transfer (Q) resolution in 1D array flattened from 2D detector
    iq_array: ndarray
        I(Q) in 1D array flattened from 2D detector
    sigma_iq_array: ndarray
        sigma I(Q) in 1D array flattened from 2D detector
    bin_centers: numpy.ndarray
        bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
    bin_edges: numpy.ndarray
        bin edges
    Returns
    -------
    ~drtsans.dataobjects.IQmod
        IQmod is a class for holding 1D binned data.

    """
    # check input
    assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

    # calculate 1/sigma^2 for multiple uses
    invert_sigma2_array = 1. / (sigma_iq_array ** 2)

    # Counts per bin: I_{k, raw} = \sum \frac{I(i, j)}{(\sigma I(i, j))^2}
    i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array * invert_sigma2_array)

    # Weight per bin: w_k = \sum \frac{1}{\sqrt{I(i, j)^2}
    w_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=invert_sigma2_array)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{w_k}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_array / w_array
    sigma_final_array = 1 / np.sqrt(w_array)

    # Calculate Q resolution of binned
    binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
    bin_q_resolution = binned_dq / i_raw_array

    # Get the final result by constructing an IQmod object defined in ~drtsans.dataobjects.
    # IQmod is a class for holding 1D binned data.
    return IQmod(intensity=i_final_array, error=sigma_final_array,
                 mod_q=bin_centers, delta_mod_q=bin_q_resolution)


def bin_intensity_into_q2d(i_of_q, qx_bins, qy_bins, method=BinningMethod.NOWEIGHT):
    """Bin I(Qx, Qy) into to new (Qx, Qy) bins

    Note: for binning parameters:
    - 'min': float or None.  If None, set to default as min(Qx) (or Qy)
    - 'max': float or None.  If None, set to default as max(Qx) (or Qy)
    - 'bins': integer as number of bins

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        class IQazimuthal(namedtuple('IQazimuthal', 'intensity error qx qy delta_qx delta_qy wavelength'))
    qx_bins : Bins
        namedtuple for arbitrary bin edges and bin centers for Qx
    qy_bins : Bins
        namedtuple for arbitrary bin edges and bin centers for Qy
    method: ~drtsans.BinningMethod
        Weighted binning or no weight binning

    Returns
    -------
    ~drtsans.dataobjects.IQazimuthal
        binned IQazimuthal

    """
    if method == BinningMethod.NOWEIGHT:
        # Calculate no-weight binning
        binned_arrays = _do_2d_no_weight_binning(i_of_q.qx, i_of_q.delta_qx, i_of_q.qy, i_of_q.delta_qy,
                                                 i_of_q.intensity, i_of_q.error, qx_bins.edges, qy_bins.edges)
    else:
        # Calculate weighed binning
        # Check input I(Q) whether it meets assumptions
        check_iq_for_binning(i_of_q)
        binned_arrays = _do_2d_weighted_binning(i_of_q.qx, i_of_q.delta_qx, i_of_q.qy, i_of_q.delta_qy,
                                                i_of_q.intensity, i_of_q.error, qx_bins.edges, qy_bins.edges)
    # END-IF-ELSE

    # construct return
    binned_intensities, binned_sigmas, binned_dqx, binned_dqy = binned_arrays

    return IQazimuthal(intensity=binned_intensities, error=binned_sigmas, qx=qx_bins.centers,
                       delta_qx=binned_dqx, qy=qy_bins.centers, delta_qy=binned_dqy)


def _do_2d_no_weight_binning(qx_array, dqx_array, qy_array, dqy_array, iq_array, sigma_iq_array,
                             qx_bin_edges, qy_bin_edges):
    """Perform 2D no-weight binning on I(Qx, Qy)

    General description of the algorithm:

      I_{i, j} = sum^{(i, j)}_k I_{k} / N_{i, j}
      sigma I_{i, j} = sqrt(sum^{(i, j)}_k sigma I_k^2) / N_{i, j}

    Parameters
    ----------
    qx_array: ndarray
        Qx array
    dqx_array: ndarray
        Qx resolution
    qy_array : ndarray
        Qy array
    dqy_array: ndarray
        Qy resolution
    iq_array: ndarray
        intensities
    sigma_iq_array: ndarray
        intensities error
    qx_bin_edges: ndarray
    qy_bin_edges

    Returns
    -------
    ndarray, ndarray, ndarray, ndarray
        intensities (n x m), sigma intensities (n x m), Qx resolution (n x m), Qy resolution (n x m)

    """
    # Number of I(q) in each target Q bin
    num_pt_array, dummy_bin_x, dummy_bin_y = np.histogram2d(qx_array, qy_array, bins=(qx_bin_edges, qy_bin_edges))

    # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
    i_raw_array, dummy_bin_x, dummy_bin_y = np.histogram2d(qx_array, qy_array, bins=(qx_bin_edges, qy_bin_edges),
                                                           weights=iq_array)

    # Square of summed uncertainties for each bin
    sigma_sqr_array, dummy_bin_x, dummy_bin_y = np.histogram2d(qx_array, qy_array, bins=(qx_bin_edges, qy_bin_edges),
                                                               weights=sigma_iq_array ** 2)

    # Q resolution: simple average
    dqx_raw_array, dummy_bin_x, dummy_bin_y = np.histogram2d(qx_array, qy_array, bins=(qx_bin_edges, qy_bin_edges),
                                                             weights=dqx_array)
    dqy_raw_array, dummy_bin_x, dummy_bin_y = np.histogram2d(qx_array, qy_array, bins=(qx_bin_edges, qy_bin_edges),
                                                             weights=dqy_array)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{Nk}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_array / num_pt_array
    sigma_final_array = np.sqrt(sigma_sqr_array) / num_pt_array
    dqx_final_array = dqx_raw_array / num_pt_array
    dqy_final_array = dqy_raw_array / num_pt_array

    return i_final_array, sigma_final_array, dqx_final_array, dqy_final_array


def _do_2d_weighted_binning(qx_array, dqx_array, qy_array, dqy_array, iq_array, sigma_iq_array,
                            x_bin_edges, y_bin_edges):
    """Perform 2D weighted binning

    General description of algorithm:

      I^{raw}_{i, j} = sum^{(i,j)}_{k} I_{k} / sigma^2(I)_k
      weight_{i, j} = sum^{(i, j)}_k 1 / sigma^2(I)_k
      I^{weight}_{i, j} = I^{raw}_{(i, j)} / weight_{i, j}
      sigma I^{weight}_{i, j} = 1 / sqrt(weight_{i, j})

    Parameters
    ----------
    qx_array : ndarray
        qx
    dqx_array : ndarray
        Qx resolution
    qy_array: ndarray
        qy
    dqy_array: ndarray
        Qy resolution
    iq_array : ndarray
        intensities
    sigma_iq_array : ndarray
        intensity errors
    x_bin_edges : ndarray
        X bin edges
    y_bin_edges
        Y bin edges

    Returns
    -------
    ndarray, ndarray, ndarray, ndarray
        binned intensities (n x m), binned sigmas (n x m), binned Qx resolution (n x m), binned Qy resolution (n x m)

    """
    # calculate 1/sigma^2 for multiple uses
    invert_sigma2_array = 1. / (sigma_iq_array ** 2)   # 1D

    # Intensities
    i_raw_2d_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                      weights=iq_array * invert_sigma2_array)  # 2D

    # dQx and dQy
    dqx_raw_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                     weights=dqx_array * invert_sigma2_array)  # 2D
    dqy_raw_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                     weights=dqy_array * invert_sigma2_array)  # 2D

    # check bins
    assert np.allclose(dummy_x, x_bin_edges, 1E-12), 'X Bin edges does not match'
    assert np.allclose(dummy_y, y_bin_edges, 1E-12), 'Y Bin edges does not match'

    # Weight per bin: w_k = \sum \frac{1}{\sqrt{I(i, j)^2}
    w_2d_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                  weights=invert_sigma2_array)  # 2D

    assert np.allclose(x_bin_edges, dummy_x, 1E-8)
    assert np.allclose(y_bin_edges, dummy_y, 1E-8)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{w_k}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_2d_array / w_2d_array
    # sigma I(Qx, Qy)
    sigma_final_array = 1 / np.sqrt(w_2d_array)
    # Qx resolution
    dqx_final_array = dqx_raw_array / w_2d_array
    # Qy resolution
    dqy_final_array = dqy_raw_array / w_2d_array

    return i_final_array, sigma_final_array, dqx_final_array, dqy_final_array
