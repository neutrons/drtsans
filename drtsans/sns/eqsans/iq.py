from __future__ import print_function

import numpy as np

from mantid.simpleapi import (ExtractSpectra, Rebin, DeleteWorkspace)
from drtsans.iq import \
    MomentumTransfer as MomentumTransferMain
from mantid.kernel import logger
from drtsans.samplelogs import SampleLogs
from drtsans.sns.eqsans import momentum_transfer

__all__ = ['prepare_momentum_transfer', 'cal_iq', 'iqxqy']

# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer(MomentumTransferMain):
    def __init__(self,
                 input_workspace=None,
                 component_name="detector1",
                 out_ws_prefix="ws"):
        super(MomentumTransfer, self).__init__(momentum_transfer, input_workspace,
                                               component_name=component_name,
                                               out_ws_prefix=out_ws_prefix)

    def __iadd__(self, other):
        """This is an overload for `+=` operator.
        Usefull for EQSANS when calculating I(qx,qy) or I(q) from several bins.

        Parameters
        ----------
        other : MomentumTransfer class

        """
        self.qx = np.concatenate((self.qx, other.qx))
        self.qy = np.concatenate((self.qy, other.qy))
        self.dqx = np.concatenate((self.dqx, other.dqx))
        self.dqy = np.concatenate((self.dqy, other.dqy))
        self.i = np.concatenate((self.i, other.i))
        self.i_sigma = np.concatenate((self.i_sigma, other.i_sigma))

        if self.component is None:
            self.component = other.component
            self.detector_dims = (self.component.dim_x, self.component.dim_y)

        return self


###############################################################################
# API
###############################################################################


def prepare_momentum_transfer(input_workspace,
                              wavelength_binning=[0.5],
                              sample_aperture=10.0,
                              prefix=None,
                              suffix="_table"):
    """ Prepare momentum transfer calculation
    :exception AssertionError: specified wave length binning parameters size is not 1 or 3
    :exception RuntimeError: wave length binni

    :param input_workspace:
    :param wavelength_binning:
    :param sample_aperture: Sample aperture diameter in mm.  It will override beamslit4 value for this aperature
    :param prefix:
    :param suffix:
    :return:
    """
    """Generates the table workspace necessary for the binning. This table contains
    unbinned Qx Qy. It is named ``prefix + suffix``, by default:
    ``input_workspace.name() + "_iqxqy_table"``

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
        The corrected Workspace
    wavelength_binning : list, optional
        This is the binning used to calculate independent I(Qi).
        In the future this will serve for I(Q) = K*I(Qi)+b.
        use: [min, step, max] or [step]
        By default [0.5], i.e. the step 0.5.
        This is the same as mantid binning format.
    sample_aperture : float, optional
        Sample aperture diameter, by default 10.0
    prefix : string, optional
        if None uses ``input_workspace.name()``, by default None
    suffix : str, optional
        The suffix for the table workspace, by default "_iqxqy_table"

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        or a :py:obj:`tuple` (:py:obj:`~mantid.api.MatrixWorkspace`, :py:obj:`~mantid.api.MatrixWorkspace`) in case
        of frame skipping datset
    """
    # Check input
    assert len(wavelength_binning) == 1 or len(wavelength_binning) == 3, \
        "wavelength_binning must be a list of 1 or 3 elements"
    # TODO - check input workspace unit

    # Add sample logs if required...
    # William: beamslit4 may not be correct. user specified sample-aperture-diameter shall override it in resolution
    #          calculation
    if sample_aperture is not None:
        input_workspace.mutableRun().addProperty('sample-aperture-diameter',
                                                 sample_aperture, 'mm', True)
    # TODO - FIXME : disabled as USELESS: geometry.source_aperture_diameter(input_workspace)

    # Identify run is 1 frame or 2 frame (skip frame)
    sl = SampleLogs(input_workspace)
    frames = []
    if bool(sl.is_frame_skipping.value):
        # case: skip frame -> 2 frames
        logger.information("This is a frame skipping data set.")
        frame1_wavelength_min = sl.wavelength_skip_min.value
        frame1_wavelength_max = sl.wavelength_skip_max.value

        frame2_wavelength_min = sl.wavelength_lead_min.value
        frame2_wavelength_max = sl.wavelength_lead_max.value

        frames.append((frame1_wavelength_min, frame1_wavelength_max))
        frames.append((frame2_wavelength_min, frame2_wavelength_max))
    else:
        # case: single frame
        frame1_wavelength_min = sl.wavelength_min.value
        frame1_wavelength_max = sl.wavelength_max.value
        frames.append((frame1_wavelength_min, frame1_wavelength_max))

    # Sanity check
    if len(wavelength_binning) > 1 and \
            bool(sl.is_frame_skipping.value):
        error_message = "The WS is frame skipping, use only the step in the wavelength_binning"
        logger.error(error_message)
        raise RuntimeError(error_message)

    # Rebin the workspace in wave length
    if prefix is None:
        prefix = input_workspace.name()
    result_wss = list()
    for index, (wavelength_min, wavelength_max) in enumerate(frames):
        # For each frame
        # set min, step, max for wave length
        if len(wavelength_binning) == 1:
            wavelength_rebinning = [wavelength_min,
                                    wavelength_binning[0], wavelength_max]
        else:
            wavelength_rebinning = wavelength_binning

        # Rebin
        ws_rebin = Rebin(InputWorkspace=input_workspace,
                         OutputWorkspace="ws_rebin",
                         Params=wavelength_rebinning)

        if len(frames) > 1:
            this_prefix = prefix + "_frame{}".format(index+1)
        else:
            this_prefix = prefix

        # Initialize a MomentumTransfer object
        mt_sum = MomentumTransfer(out_ws_prefix=this_prefix)

        # Calculate momentum transfer and set up the I(Q)-TableWorkspace
        bins = ws_rebin.readX(0)
        for bin_index in range(len(bins)-1):
            # Calculate Q for a single bin of wave length
            ws_extracted = ExtractSpectra(InputWorkspace=ws_rebin,
                                          XMin=bins[bin_index],
                                          XMax=bins[bin_index+1])

            # TODO FIXME #214: The commented codes seem going nowhere
            # wavelength_mean = (bins[bin_index] + bins[bin_index+1]) / 2.0
            # runObj = ws_extracted.mutableRun()
            # runObj.addProperty('wavelength', float(wavelength_mean),
            #                    'Angstrom', True)
            # runObj.addProperty('wavelength-spread', 0.2, 'Angstrom', True)

            # core calculation:
            mt_extracted = MomentumTransfer(ws_extracted)
            # Append to final Momentum
            mt_sum += mt_extracted
        # END-FOR

        # Create Q2D table workspace and add to return TableWorkspace list
        result_wss.append(mt_sum.q2d(suffix=suffix))

        # Clean
        DeleteWorkspace(ws_rebin)
        DeleteWorkspace(ws_extracted)
    # END-FOR (frame)

    return result_wss


def cal_iq(input_table_workspace, bins=100, log_binning=False, suffix="_iq"):
    """
    Creates a WS named: input_workspace.name() + "_iq"

    Parameters
    ----------
    input_table_workspace : ~mantid.api.ITableWorkspace

    bins : int or sequence of scalars, optional
        See `scipy.stats.binned_statistic
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binned_statistic.html>`_.
        If ``bins`` is an int, it defines the number of equal-width bins in
        the given range (10 by default).  If ``bins`` is a sequence, it
        defines the bin edges, including the rightmost edge, allowing for
        non-uniform bin widths.  Values in ``x`` that are smaller than lowest
        bin edge areassigned to bin number 0, values beyond the highest bin
        are assigned to ``bins[-1]``.  If the bin edges are specified,
        the number of bins will be, (nx = len(bins)-1).
    log_binning : bool, optional
        if True bins must be an integer, by default False
    Returns
    -------
    ~mantid.api.MatrixWorkspace
        Only a WS with a single spectra
    """

    mt = MomentumTransfer(
        input_table_workspace,
        out_ws_prefix=input_table_workspace.name().rstrip('_table'))

    if log_binning and isinstance(bins, int):
        q = np.sqrt(np.square(mt.qx) + np.square(mt.qy))
        bins = np.logspace(np.log10(np.min(q)), np.log10(np.max(q)), num=bins)

    _, ws = mt.bin_into_q1d(bins=bins, suffix=suffix)
    return ws


def _linear_log_array(arr, n_bins):
    """Transforms a linear binning array into a log
    The input array must be sorted, for example:
    [-5, ..., -1, 0, 1, ..., 6]
    is valid.

    Parameters
    ----------
    arr : array
        array with original binning with negative values
    n_bins : int
        number of bins in the returned array

    Returns
    -------
    array
        Logarithmic binning array
    """

    negative_bins = np.array([arr[0], np.where(np.array(arr) < 0,
                                               arr, -np.inf).max()])
    positive_bins = np.array([np.where(np.array(arr) > 0,
                                       arr, np.inf).min(), arr[-1]])

    negative_bins *= -1
    negative_bins = np.flip(negative_bins, axis=0)

    bins_negative = np.geomspace(negative_bins[0], negative_bins[1], n_bins//2)
    bins_negative *= -1
    bins_negative = np.flip(bins_negative, axis=0)

    bins_positive = np.geomspace(positive_bins[0], positive_bins[1], n_bins//2)
    bins = np.hstack((bins_negative, bins_positive))

    return bins


def iqxqy(input_table_workspace, bins=100, log_binning=False, suffix='_iqxqy'):
    """
    Only used to visulize Qx vs Qy.
    Creates a WS named: ``str(input_table_workspace) + '_iqxqy'``

    Parameters
    ----------
    input_table_workspace : ~mantid.api.ITableWorkspace

    bins : int or array_like or [int, int] or [array, array], optional
        The bin specification:
            - If int, the number of bins for the two dimensions (nx=ny=bins).
            - If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
            - If [int, int], the number of bins in each dimension (nx, ny = bins).
            - If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).
            - A combination [int, array] or [array, int], where int is the number of bins and array is the bin edges.

        If :py:obj:`None` the detector dimensions are given:
        bins = [det.dim_x, det.dim_y]

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        For visualization only.
    """

    mt = MomentumTransfer(
        input_table_workspace,
        out_ws_prefix=input_table_workspace.name().rstrip('_table'))

    if log_binning:
        # number of bins
        num_x = None
        num_y = None
        if isinstance(bins, int):
            num_x = bins
            num_y = bins
        elif isinstance(bins, list) and len(bins) == 2:
            num_x = bins[0]
            num_y = bins[1]
        else:
            logger.error("Please use bins as an int or a list with two "
                         "arguments. E.g.: bins=100 or bins=[120, 80]")
            return

        if num_x is not None and num_y is not None:
            # first let's find where the negatives start and end
            bins_qx = _linear_log_array(sorted(mt.qx), num_x)
            bins_qy = _linear_log_array(sorted(mt.qy), num_y)
            bins = [bins_qx, bins_qy]

    _, ws = mt.bin_into_q2d(bins=bins, suffix=suffix)

    return ws


def iq_wedge(input_table_workspace,  phi_0=0,
             phi_aperture=30, bins=100, statistic='mean',
             suffix="_wedge_iq"):
    """
    Wedge calculation and integration

    Calculates: I(Q) and Dq
    The ws_* input parameters are the output workspaces from bin_into_q2d

    Parameters
    ----------
    input_table_workspace : TableWorkspace

    phi_0 : int, optional
        Where to start the wedge, by default 0
    phi_aperture : int, optional
        Aperture of the wedge, by default 30
    bins : int or sequence of scalars, optional
        See `scipy.stats.binned_statistic`.
        If `bins` is an int, it defines the number of equal-width bins in
        the given range (10 by default).  If `bins` is a sequence, it
        defines the bin edges, including the rightmost edge, allowing for
        non-uniform bin widths.  Values in `x` that are smaller than lowest
        bin edge areassigned to bin number 0, values beyond the highest bin
        are assigned to ``bins[-1]``.  If the bin edges are specified,
        the number of bins will be, (nx = len(bins)-1).
    statistic : str, optional
        See `scipy.stats.binned_statistic`.
        The statistic to compute, by default 'mean'
        The following statistics are available:
        * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
        * 'std' : compute the standard deviation within each bin. This
            is implicitly calculated with ddof=0.
        * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
        * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
        * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
        * 'min' : compute the minimum of values for points within each bin.
            Empty bins will be represented by NaN.
        * 'max' : compute the maximum of values for point within each bin.
            Empty bins will be represented by NaN.
        * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.
    out_ws_prefix : str, optional
        The prefix of the workspace created in Mantid, by default "ws"

    Returns
    -------
        Workspace2D with one spectrum: I(q)

    """

    mt = MomentumTransfer(
        input_table_workspace,
        out_ws_prefix=input_table_workspace.name().rstrip('_table'))

    _, ws = mt.bin_wedge_into_q1d(phi_0=phi_0, phi_aperture=phi_aperture,
                                  bins=bins, statistic=statistic,
                                  suffix=suffix)
    return ws


def iq_annular(input_table_workspace, q_min=0.001, q_max=0.4,
               bins=100, statistic='mean', suffix="_annular_iq"):
    """
    Wedge calculation and integration

    Calculates: I(Q) and Dq
    The ws_* input parameters are the output workspaces from bin_into_q2d

    Parameters
    ----------
    input_table_workspace : TableWorkspace

    q_min : float, optional
        , by default
    q_max : float, optional
        , by default
        Aperture of the wedge, by default 30
    bins : int or sequence of scalars, optional
        See `scipy.stats.binned_statistic`.
        If `bins` is an int, it defines the number of equal-width bins in
        the given range (10 by default).  If `bins` is a sequence, it
        defines the bin edges, including the rightmost edge, allowing for
        non-uniform bin widths.  Values in `x` that are smaller than lowest
        bin edge areassigned to bin number 0, values beyond the highest bin
        are assigned to ``bins[-1]``.  If the bin edges are specified,
        the number of bins will be, (nx = len(bins)-1).
    statistic : str, optional
        See `scipy.stats.binned_statistic`.
        The statistic to compute, by default 'mean'
        The following statistics are available:
        * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
        * 'std' : compute the standard deviation within each bin. This
            is implicitly calculated with ddof=0.
        * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
        * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
        * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
        * 'min' : compute the minimum of values for points within each bin.
            Empty bins will be represented by NaN.
        * 'max' : compute the maximum of values for point within each bin.
            Empty bins will be represented by NaN.
        * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.
    out_ws_prefix : str, optional
        The prefix of the workspace created in Mantid, by default "ws"

    Returns
    -------
        Workspace2D with one spectrum: I(q)

    """

    mt = MomentumTransfer(
        input_table_workspace,
        out_ws_prefix=input_table_workspace.name().rstrip('_table'))

    _, ws = mt.bin_annular_into_q1d(q_min=q_min, q_max=q_max, bins=bins,
                                    statistic=statistic, suffix=suffix)
    return ws
