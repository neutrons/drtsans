# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/dataobjects.py
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/docs/drtsans/dataobjects.rst
from drtsans.dataobjects import IQazimuthal, IQmod
from enum import Enum
import numpy as np
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/determine_bins.py
from drtsans.determine_bins import determine_1d_linear_bins, determine_1d_log_bins, BinningParams
# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')

__all__ = ['bin_intensity_into_q1d', 'select_i_of_q_by_wedge',
           'bin_annular_into_q1d', 'bin_intensity_into_q2d', 'BinningMethod', 'check_iq_for_binning',
           'determine_1d_linear_bins', 'determine_1d_log_bins', 'BinningParams']


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

    # Check input I(Q) whether it meets assumptions
    check_iq_for_binning(i_of_q)

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
        binned_iq = _do_1d_weighted_binning(theta_array[allowed_q_index],
                                            dq_array[allowed_q_index],
                                            i_of_q.intensity[allowed_q_index],
                                            i_of_q.error[allowed_q_index],
                                            theta_bins.centers, theta_bins.edges)
    else:
        # not supported case
        raise RuntimeError('Binning method {} is not recognized'.format(method))

    return binned_iq


def _do_1d_no_weight_binning(q_array, dq_array, iq_array, sigmaq_array, bin_centers, bin_edges):
    """ Bin I(Q) by given bin edges and do no-weight binning

    This method implements equation 11.34, 11.35 and 11.36 in master document.

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

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

    This method implements equation 11.22, 11.23 and 11.24 in master document for 1-dimensional Q

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

    General description of algorithm:

    Equation 11.22
    I(Q') = sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2) /
            sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)

    Equation 11.23
    sigmaI(Q') = sqrt(sum_{Q, lambda}^{K} (sigma(Q, lambda / sigma(Q, lambda)^2)^2) /
                 sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
               = sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)) /
                 sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
               = 1 / sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2))

    Equation 11.24
    sigmaQ(Q') = sum_{Q, lambda}^{K}(sigmaQ(Q, lambda)/sigma^2(Q, lambda)^2) /
                 sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)

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
    # Check input
    assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

    # Calculate 1/sigma^2 for multiple uses
    invert_sigma2_array = 1. / (sigma_iq_array ** 2)

    # Histogram on 1/sigma^2, i.e., nominator part in Equation 11.22, 11.23 and 11.24
    # sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
    w_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=invert_sigma2_array)

    # Calculate Equation 11.22: I(Q)
    #  I(Q') = sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2) /
    #              sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
    # denominator in Equation 11.22: sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2)
    i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array * invert_sigma2_array)
    # denominator divided by nominator (11.22)
    i_final_array = i_raw_array / w_array

    # Calculate equation 11.23: sigmaI(Q)
    # sigmaI(Q') = sqrt(sum_{Q, lambda}^{K} (sigma(Q, lambda / sigma(Q, lambda)^2)^2) /
    #              sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
    #            = sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)) /
    #              sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
    #             = 1 / sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2))
    # Thus histogrammed sigmaI can be obtained from histogrammed invert_sigma2_array directly
    sigma_final_array = 1 / np.sqrt(w_array)

    # Calculate equation 11.24:  sigmaQ (i.e., Q resolution)
    # sigmaQ(Q') = sum_{Q, lambda}^{K}(sigmaQ(Q, lambda)/sigma^2(Q, lambda)^2) /
    #              sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
    # denominator in Equation 11.24: sum_{Q, lambda}^{K}(sigmaQ(Q, lambda)/sigma^2(Q, lambda)^2)
    binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
    # denominator divided by nominator (11.24)
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
    # Check input I(Q) whether it meets assumptions
    check_iq_for_binning(i_of_q)

    if method == BinningMethod.NOWEIGHT:
        # Calculate no-weight binning
        binned_arrays = _do_2d_no_weight_binning(i_of_q.qx, i_of_q.delta_qx, i_of_q.qy, i_of_q.delta_qy,
                                                 i_of_q.intensity, i_of_q.error, qx_bins.edges, qy_bins.edges)
    else:
        # Calculate weighed binning
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

    I(x', y')      = sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2) /
                     sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
    sigmaI(x', y') = sqrt(sum_{x, y, lambda}^{K} (sigma(x, y, lambda / sigma(x, y, lambda)^2)^2) /
                     sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
                   = sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)) /
                     sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
                   = 1 / sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2))
    sigmaQ(x', y') = sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2) /
                     sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)

    where K is the set of (x, y, sigma) such that (x, y, sigma) in the same Q_bin

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
    # Calculate 1/sigma^2 for multiple uses
    invert_sigma2_array = 1. / (sigma_iq_array ** 2)   # 1D

    # Histogram on 1/sigma^2, i.e., nominator part in Equation 11.22, 11.23 and 11.24
    # sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
    w_2d_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                  weights=invert_sigma2_array)  # 2D

    # Calculate Equation 11.22: I(Qx, Qy)
    # I(x', y') = sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2) /
    #             sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
    # denominator in Equation 11.22: sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2)
    i_raw_2d_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                      weights=iq_array * invert_sigma2_array)  # 2D
    # denominator divided by nominator (11.22)
    i_final_array = i_raw_2d_array / w_2d_array

    # Calculate equation 11.23: sigmaI(Q)
    # sigmaI(x', y') = sqrt(sum_{x, y, lambda}^{K} (sigma(x, y, lambda / sigma(x, y, lambda)^2)^2) /
    #                  sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
    #                = sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)) /
    #                sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
    #                = 1 / sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2))
    # Thus histogrammed sigmaI can be obtained from histogrammed invert_sigma2_array directly
    sigma_final_array = 1 / np.sqrt(w_2d_array)

    # Calculate equation 11.24:  sigmaQx and sigmaQy (i.e., Q resolution)
    # sigmaQ(x', y') = sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2) /
    #                  sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
    # denominator in Equation 11.24: sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2)
    dqx_raw_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                     weights=dqx_array * invert_sigma2_array)  # 2D
    dqy_raw_array, dummy_x, dummy_y = np.histogram2d(qx_array, qy_array, bins=(x_bin_edges, y_bin_edges),
                                                     weights=dqy_array * invert_sigma2_array)  # 2D
    # denominator divided by nominator (11.24)
    dqx_final_array = dqx_raw_array / w_2d_array  # dQx
    dqy_final_array = dqy_raw_array / w_2d_array  # dQy

    return i_final_array, sigma_final_array, dqx_final_array, dqy_final_array
