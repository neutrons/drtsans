from string import Template

import numpy as np
from scipy import stats

import mantid
from mantid.simpleapi import CreateEmptyTableWorkspace, CreateWorkspace
from mantid.api import AnalysisDataService
from drtsans.momentum_transfer_factory import calculate_q_dq
from drtsans.detector import Component
import collections
from enum import Enum

# from mantid.kernel import logger
# from ast import literal_eval
# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')

"""Proposed API
- binning Q to 1D:               q_vector, dq_vector = bin_into_q1d(wl_ws, bins, statistic)
- binning Q to 2D:               qx_vector, dqx_vector, qy_vector, dqy_vector = bin_into_q2d(wl_ws, bins, statistic)
- binning Q to 1D with annular:  q_vector, dq_vector = bin_annular_into_q1d(wl_ws, bins, statistic)
- binning Q to 1D with wedge:    q_vector, dq_vector = bin_wedge_into_q1d(wl_ws, bins, statistic)
- exporting result to table:     tables = export_to_tables(q_vectors, dq_vectors)
"""

# Define structure for I(Q) for Q, dQ, I(Q), Sigma_I(Q)
IofQ = collections.namedtuple('IofQ', 'q dq i sigma')

# Define structure (namedtuple) for binning parameters: min, max, number of bins
# bins shall be integer as number of bins
BinningParams = collections.namedtuple('BinningParams', 'min max bins')


class BinningMethod(Enum):
    """
    Binning method
    """
    NOWEIGHT = 1
    WEIGHTED = 2


def bin_iq_into_linear_q1d(wl_ws, bins, q_min=0, q_max=None, instrument=None):
    """
    Binning I(Q) from a mono
    Parameters
    ----------
    wl_ws :  ~mantid.api.MatrixWorkspace
        MatrixWorkspace in wave length
    bins: integer
        number of bins
    q_min: float
        minimum Q left bin edge value
    q_max: float or None
        maximum Q right bin edge value. None as the default such that Qmax will be determined from workspace
            automatically
    instrument: String or None
        Instrument type (tof, mono) only used for testing purpose. For real SANS instrument, it will be determined
        automatically from instrument name
    Returns
    -------
    IofQ
        named tuple for Q, dQ, I(Q), sigma_I(Q)
        numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
        Q, dQ, I, dI
        Q, Q resolution, I, uncertainty of I
    """
    # Initialize the calculator
    calculator = IofQCalculator(wl_ws, instrument_type=instrument)

    # retrieve Qmax
    if q_max is None:
        q_max = calculator.get_q_max()

    # bin
    binned_i_q = calculator.bin_into_linear_q1d(bins, q_min, q_max)

    return binned_i_q


def bin_iq_into_logarithm_q1d(wl_ws, bins_per_decade, q_min=0.001, q_max=1.0, instrument=None):
    """
    Binning I(Q) with a logarithm bin
    Parameters
    ----------
    wl_ws :  ~mantid.api.MatrixWorkspace
        MatrixWorkspace in wave length
    bins_per_decade: integer
        number of bins
    q_min: float
        minimum Q left bin edge value
    q_max: float or None
        maximum Q right bin edge value. None as the default such that Qmax will be determined from workspace
            automatically
    instrument: String or None
        Instrument type (tof, mono) only used for testing purpose. For real SANS instrument, it will be determined
        automatically from instrument name
    Returns
    -------
    IofQ
        named tuple for Q, dQ, I(Q), sigma_I(Q)
        numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
        Q, dQ, I, dI
        Q, Q resolution, I, uncertainty of I
    """
    # Initialize the calculator
    calculator = IofQCalculator(wl_ws, instrument_type=instrument)

    # bin
    binned_i_q = calculator.bin_into_log_q1d(bins_per_decade, q_min, q_max)

    return binned_i_q


def bin_into_q2d(wl_ws, bins, suffix):
    """
    :param wl_ws: List of workspaces (names) in binned wave length space
    :param bins: Iterable for range and bin size of Qx and Qy
    :param suffix: suffix for output workspace
    :return:
    """
    calculator = IofQCalculator(wl_ws)

    return calculator.bin_into_q2d(bins=bins, suffix=suffix)


def bin_iq_into_linear_q2d(i_q, qx_bin_params, qy_bin_params, method=BinningMethod.NOWEIGHT):
    """Bin I(Qx, Qy) into to new (Qx, Qy) bins

    Note: for binning parameters:
    - 'min': float or None.  If None, set to default as min(Qx) (or Qy)
    - 'max': float or None.  If None, set to default as max(Qx) (or Qy)
    - 'bins': integer as number of bins

    Parameters
    ----------
    i_q: namedtuple
        "i": intensity, "qx": qx, "qy": qy, "dqx": dqx, "dqy", dqy
    qx_bin_params: BinningParams
        binning parameters for Qx
    qy_bin_params: BinningParams
        binning parameters for Qy
    method: BinningMethod
        Weighted binning or no weight binning

    Returns
    -------

    """
    # Calculate Qx and Qy bin size
    qx_bin_size = _determine_linear_bin_size(i_q.qx, qx_bin_params.min, qx_bin_params.bins, qx_bin_params.max)
    qy_bin_size = _determine_linear_bin_size(i_q.qy, qy_bin_params.min, qy_bin_params.bins, qy_bin_params.max)
    print(qx_bin_size, qy_bin_size)

    # Calculate histogram

    return


def _determine_linear_bin_size(x_array, min_x, num_bins, max_x):
    """Determine linear bin size

    This is adopted by bin I(Qx, Qy)

    Parameters
    ----------
    x_array: ndarray
        Value X
    min_x: float
        minimum X. None as default x_array.min()
    num_bins: integer
        number of bins
    max_x: float
        maximum X. None as default x_array.max()

    Returns
    -------

    """
    # Determine min X and max X
    if min_x is None:
        min_x = np.min(x_array)
    if max_x is None:
        max_x = np.max(x_array)

    # Calculate delta
    if num_bins <= 1:
        raise RuntimeError('Number of bins cannot be less than 2')

    delta_x = (max_x - min_x) / (num_bins - 1.)

    return delta_x


def bin_wedge_into_q1d(wl_ws, phi_0=0, phi_aperture=30, bins=100,
                       statistic='mean', suffix="_wedge_iq"):
    """
    Wedge calculation and integration

        Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------
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
        suffix : str, optional
            The prefix of the workspace created in Mantid, by default "ws"

        Returns
        -------
        workspaces list
    """
    calculator = IofQCalculator(wl_ws)

    return calculator.bin_wedge_into_q1d(phi_0, phi_aperture, bins, statistic, suffix)


def bin_annular_into_q1d(wl_ws,
                         q_min=0.001,
                         q_max=0.4,
                         bins=100,
                         statistic='mean',
                         suffix="_annular_iq"):
    """
    Wedge calculation and integration

    Calculates: I(Q) and Dq
    The ws_* input parameters are the output workspaces from bin_into_q2d

    Parameters
    ----------
    wl_ws : list of workspaces (reference or string as name)
    q_min : float, optional
        , by default
    q_max : float, optional
        , by default
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
    suffix : str, optional
        The prefix of the workspace created in Mantid, by default "ws"

    Returns
    -------
    tuple
        (workspace name, workspace)

    """
    calculator = IofQCalculator(wl_ws)

    return calculator.bin_annular_into_q1d(q_min, q_max, bins, statistic, suffix)


def export_i_q_to_table(i_of_q, table_ws_name, detector_dims, DETECTOR_DIMENSIONS_TEMPLATE):
    """
    Export binned I(Q) to table (workspace)
    Returns
    -------

    """
    def _create_table_ws(table_ws_name, detector_dim):
        """
        Create a Mantid TableWorkspace containing raw Qx, Qy, dQx, dQy, I(Q) and Sigma(Q)
        Parameters
        ----------
        prefix: String
            prefix of the output TableWorkspace
        suffix: String
            suffix of the output TableWorkspace

        Returns
        -------
        mantid.api.ITableWorkspace
            TableWorkspace containing Qx, Qy, dQx, dQy, I and sigma(I)s
        """
        # Create empty table
        table_iq = CreateEmptyTableWorkspace(OutputWorkspace=table_ws_name)
        # Add columns for Q (2D)
        table_iq.addColumn(type="float", name="Qx")
        table_iq.addColumn(type="float", name="Qy")
        table_iq.addColumn(type="float", name="dQx")
        table_iq.addColumn(type="float", name="dQy")
        table_iq.addColumn(type="float", name="I")
        table_iq.addColumn(type="float", name="Sigma(I)")

        return table_iq

    # Create workspace
    table_iq = _create_table_ws(table_ws_name)

    # Add each (Qx, Qy, I(Qx, Qy)) to table workspace
    for qx_i, qy_i, dqx_i, dqy_i, i_i, i_sigma_i in zip(i_of_q.qx.tolist(), i_of_q.qy.tolist(),
                                                        i_of_q.dqx.tolist(), i_of_q.dqy.tolist(),
                                                        i_of_q.i_q.tolist(), i_of_q.sigma_i_q.tolist()):
        new_row = {'Qx': qx_i,
                   'Qy': qy_i,
                   'dQx': dqx_i,
                   'dQy': dqy_i,
                   "I": i_i,
                   "Sigma(I)": i_sigma_i
                   }
        table_iq.addRow(new_row)
    # END-FOR

    # Add comment
    template = Template(DETECTOR_DIMENSIONS_TEMPLATE)
    table_iq.setComment(
            template.substitute(dim_x=detector_dims[0],
                                dim_y=detector_dims[1]))

    return table_iq
    # data = input_workspace.toDict()
    # self.qx = np.array(data['Qx'])
    # self.qy = np.array(data['Qy'])
    # self.dqx = np.array(data['dQx'])
    # self.dqy = np.array(data['dQy'])
    # self.i = np.array(data["I"])
    # self.i_sigma = np.array(data["Sigma(I)"])
    #
    # # Gets the detector dimensions: "detector_dimensions=($dim_x,$dim_y)"
    # if self.DETECTOR_DIMENSIONS_TEMPLATE.split("=")[0] != \
    #         input_workspace.getComment().split("=")[0]:
    #     logger.error("Can not get the detector dimensions!")
    # else:
    #     self.detector_dims = literal_eval(
    #         input_workspace.getComment().split("=")[1])


class IofQCalculator(object):
    """
    Momentum Transfer class
    Olds arrays for qx, qy, dqx, dqy, i, i_sigma
    q1d and q2d are calculated from there arrays
    """

    # For now the detector dims in Table Workspace will be in the comment
    # as below (there is no other form of adding key=value to table WS)
    DETECTOR_DIMENSIONS_TEMPLATE = "detector_dimensions=($dim_x,$dim_y)"

    qx, qy, dqx, dqy, i, i_sigma = np.empty(0), np.empty(0), np.empty(0), \
        np.empty(0), np.empty(0), np.empty(0)
    prefix = None
    component = None
    detector_dims = None

    def __init__(self, input_workspace, component_name="detector1",
                 out_ws_prefix="ws", instrument_type=None):
        """
        Initialization of I(Q) calculator
        Parameters
        ----------
        input_workspace: String or  ~mantid.api.MatrixWorkspace
            Workspace name or workspace instance
        component_name: String
            component name for detector information coming from
        out_ws_prefix: String
            prefix for output workspace
        instrument_type: String
            instrument type: mono or tof. used for testing purpose
        """
        # Set up class variable
        # input wave length workspace
        if not isinstance(input_workspace, mantid.api.MatrixWorkspace):
            raise RuntimeError('Input workspace {} must be a MatrixWorkspace'.format(input_workspace))
        self._wl_ws = AnalysisDataService.retrieve(str(input_workspace))
        # output workspace's prefix
        self._prefix = out_ws_prefix

        # get detector dimension
        component = Component(input_workspace, component_name)
        self._detector_dims = component.dim_x, component.dim_y

        # Calculate q, qx, qy, dqx, dqy: N x M array where N is number of histograms in input workspace
        # It is not necessary to consider monitor or masked detectors (pixels) for Q and dQ
        self._q_dq = calculate_q_dq(input_workspace, instrument_type=instrument_type)

        # Extract I(Q) from input workspace: N x M array
        # such that self._i_q[n, m] corresponds to self._q_dq[n, m]
        self._i_q = input_workspace.extractY().ravel()
        self._i_q_sigma = input_workspace.extractE().ravel()

        # Mask monitors and detectors pixels
        masked_pixels = component.masked_ws_indices()
        monitor_pixels = component.monitor_indices()  # TODO FIXME - check monitor_indices implemented?
        self._mask_pixels(masked_pixels, monitor_pixels)\

        # 1D arrays
        # self.qx, self.qy, self.dqx, self.dqy = \
        #     self.resolution.calculate_q_dq(input_workspace)
        # # Ravel just in case! For EQSANS for example is needed!
        # self.qx, self.qy, self.dqx, self.dqy = \
        #     self.qx.ravel(), self.qy.ravel(), self.dqx.ravel(), \
        #     self.dqy.ravel()
        # # Get rid of the monitors; from 49154 to 49152 spectra
        # self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
        #     self._remove_monitors(d) for d in
        #     [self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]
        # ]
        #
        # # Create numpy mask arrays with the masked pixels
        # self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
        #     self._mask_pixels(d, masked_pixels) for d in
        #     [self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]
        # ]

        return

    def get_q_max(self):
        """
        Get Q max
        Returns
        -------
        float
            Q max
        """
        return np.max(self._q_dq.q)

    def _mask_pixels(self, masked_pixels, monitor_pixels):
        """
        Mask pixels
        Parameters
        ----------
        data
        masked_pixels

        Returns
        -------

        """
        # Mask I(Q)
        try:
            self._i_q = np.ma.MaskedArray(self._i_q,
                                          masked_pixels,
                                          dtype=np.float,
                                          fill_value=np.nan)

            # Mask sigma I(Q)
            self._i_q_sigma = np.ma.MaskedArray(self._i_q_sigma,
                                                masked_pixels,
                                                dtype=np.float,
                                                fill_value=np.nan)
        except np.ma.core.MaskError:
            # TODO FIXME - Masking is not correct!
            pass

        return

    def bin_into_q2d(self, bins=None, suffix="_iqxqy"):
        """Bin the data into Q 2D for visualization only!
        No uncertainties data!!!

        Parameters
        ----------

        bins : int or array_like or [int, int] or [array, array], optional
            The bin specification:
                - If int, the number of bins for the two dimensions
                (nx=ny=bins).
                - If array_like, the bin edges for the two dimensions
                (x_edges=y_edges=bins).
                - If [int, int], the number of bins in each dimension
                (nx, ny = bins).
                - If [array, array], the bin edges in each dimension
                (x_edges, y_edges = bins).
                - A combination [int, array] or [array, int], where int is the
                number of bins and array is the bin edges.
            If None the detector dimensions are given:
            bins = [det.dim_x, det.dim_y]

        Returns
        -------
        string, Workspace2D
            workspace name, workspace reference

        """
        # Set up default bins:
        # TODO - confirm component.dim_x and component.dim_y is the Qx, Qy binning range
        if bins is None:
            bins = self._detector_dims

        # TODO - Need to verify the Nan (masked) I(Q) can do its job
        # Number of bins in Qx Qy is the number of pixels in X and Y
        # Bin for number of pixels in each bin
        counts_qx_qy_weights, qx_bin_edges, qy_bin_edges = np.histogram2d(self._q_dq.qx, self._q_dq.qy, bins=bins)
        # Bin for I(Q)
        counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(self.qx, self.qy, bins=bins, weights=self.i)
        # Normalize
        counts_qx_qy /= counts_qx_qy_weights
        # Convert from bin edgets to bin centers
        qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0

        # When doing histogram2d the masks are gone
        iqxqy_ws = CreateWorkspace(DataX=np.tile(qx_bin_edges,
                                                 len(qy_bin_centers)),
                                   DataY=np.fliplr(counts_qx_qy.T),
                                   NSpec=len(qy_bin_centers),
                                   UnitX='MomentumTransfer',
                                   VerticalAxisUnit='MomentumTransfer',
                                   VerticalAxisValues=qy_bin_centers,
                                   OutputWorkspace=self.prefix + suffix)

        return iqxqy_ws.name(), iqxqy_ws

    # TODO - need to verify whether NaN can work for masked pixels
    @staticmethod
    def _bin_intensity_into_q1d(q, dq, i, i_sigma, prefix, bins=100,
                                statistic='mean', suffix="_iq"):
        """
        Bin I(Q) by scaler Q (1D) with given number of bins
        Parameters
        ----------
        q
        dq
        i
        i_sigma
        prefix
        bins
        statistic
        suffix

        Returns
        -------

        """
        # Bin I(Q)
        intensity_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, i,
            statistic=statistic, bins=bins)

        # Bin Sigma(I)
        sigma_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, i_sigma,
            statistic=lambda array_1d: np.sqrt(np.sum(np.square(array_1d))
                                               ) / len(array_1d),
            bins=bins)

        # Bin for weight (pixels in each bin)
        _, dq_bin_edges, _ = \
            stats.binned_statistic(dq, i, statistic=statistic, bins=bins)

        # d(Q) for bin center
        dq_bin_centers = (dq_bin_edges[1:] + dq_bin_edges[:-1]) / 2.0

        # Create MomentumTransfer Workspace2D
        iq = CreateWorkspace(
            DataX=np.array([q_bin_edges]),
            DataY=np.array([intensity_statistic]),
            DataE=np.array([sigma_statistic]),
            Dx=dq_bin_centers,  # bin centers!!
            NSpec=1,
            UnitX='MomentumTransfer',
            YUnitLabel='Counts',
            OutputWorkspace=prefix + suffix)

        return iq.name(), iq

    def bin_into_linear_q1d(self, bins, q_min, q_max):
        """
        Bin I(Q) by linear binning with linear binning
        Parameters
        ----------
        bins : integer
            number of bins
        q_min: float
            minimum Q (center) value of Q bins
        q_max: float or None
            maximum Q (center) value of Q bins.  None as the default as max Q

        Returns
        -------
        IofQ
            named tuple for Q, dQ, I(Q), sigma_I(Q)
        """
        # Determine bin edges
        # delta_q = (q_max - q_min) / bins
        # bin_edges = np.arange(bins + 1).astype('float') * delta_q + q_min
        # bin_centers = (bin_edges[1:] - bin_edges[:-1]) * 0.5
        bin_centers, bin_edges = self.determine_linear_bin_edges(q_min, q_max, bins)

        # Bin
        binned_i_q = self.weighted_binning(self._q_dq.q, self._q_dq.dq, self._i_q, bin_centers, bin_edges)

        return binned_i_q

    def bin_into_log_q1d(self, bins_per_decade, q_min, q_max):
        """

        Parameters
        ----------
        bins_per_decade
        q_min
        q_max

        Returns
        -------

        """
        bin_centers, bin_edges = self.determine_log_bin_edges(q_min, q_max, bins_per_decade)

        # Bin
        binned_i_q = self.weighted_binning(self._q_dq.q, self._q_dq.dq, self._i_q, bin_centers, bin_edges)

        return binned_i_q

    def bin_wedge_into_q1d(self, phi_0=0, phi_aperture=30, bins=100,
                           statistic='mean', suffix="_wedge_iq"):
        """
        Wedge calculation and integration

        Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------
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
        tuple
            tuple (workspace name, workspace)

        """

        angle_rad = np.arctan2(self.qy, self.qx)

        # plt.scatter(mt.qx, mt.qy, c=angle_rad)
        # plt.colorbar()
        # plt.show()

        # # This is just to show the angle
        # CreateWorkspace(
        #     DataX=self.qx,
        #     DataY=angle_grid,
        #     DataE=np.sqrt(angle_grid),
        #     NSpec=256,
        #     UnitX='MomentumTransfer',
        #     VerticalAxisUnit='MomentumTransfer',
        #     VerticalAxisValues=self.qy,
        #     OutputWorkspace=out_ws_prefix+"_angle",
        # )

        # Let's work in radians
        phi_0_rad = np.deg2rad(phi_0)
        phi_aperture_rad = np.deg2rad(phi_aperture)

        phi_aperture_min_rad = phi_0_rad - phi_aperture_rad / 2
        phi_aperture_max_rad = phi_0_rad + phi_aperture_rad / 2
        # opposite 180 degrees apart
        phi_aperture_min_pi_rad = phi_aperture_min_rad + np.pi
        phi_aperture_max_pi_rad = phi_aperture_max_rad + np.pi

        condition1 = (angle_rad > phi_aperture_min_rad) & \
            (angle_rad < phi_aperture_max_rad)
        # make angle > np.pi varying between np.pi and 2*np.pi, rather than the
        # initial -np.pi to np.pi
        angle_rad[angle_rad < 0] = 2 * np.pi + angle_rad[angle_rad < 0]
        condition2 = (angle_rad > phi_aperture_min_pi_rad) & \
            (angle_rad < phi_aperture_max_pi_rad)
        # True where the wedge is located, otherwise false
        condition = condition1 | condition2

        # # This is just to show the condition
        # CreateWorkspace(
        #     DataX=qx_bin_edges_grid,
        #     DataY=condition,
        #     DataE=np.sqrt(angle_grid),
        #     NSpec=256,
        #     UnitX='MomentumTransfer',
        #     VerticalAxisUnit='MomentumTransfer',
        #     VerticalAxisValues=qy_bin_centers,
        #     OutputWorkspace=out_ws_prefix+"_condition",
        # )

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        q = q[condition]

        i = self.i[condition]
        i_sigma = self.i_sigma[condition]

        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))
        dq = dq[condition]

        return IofQCalculator._bin_intensity_into_q1d(q, dq, i, i_sigma,
                                                      self.prefix, bins,
                                                      statistic, suffix)

    def bin_annular_into_q1d(self,
                             q_min=0.001,
                             q_max=0.4,
                             bins=100,
                             statistic='mean',
                             suffix="_annular_iq"):
        """
        Wedge calculation and integration

        Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------
        q_min : float, optional
            , by default
        q_max : float, optional
            , by default
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
        tuple
            (workspace name, workspace)

        """

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        condition = (q >= q_min) & (q <= q_max)

        q = q[condition]
        i = self.i[condition]
        i_sigma = self.i_sigma[condition]

        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))
        dq = dq[condition]

        return IofQCalculator._bin_intensity_into_q1d(q, dq, i, i_sigma,
                                                      self.prefix, bins,
                                                      statistic, suffix)

    @staticmethod
    def weighted_binning(q_array, dq_array, iq_array, sigmaq_array, bin_centers, bin_edges):
        """ Bin I(Q) by given bin edges and do weighted binning
        Parameters
        ----------
        q_array: ndarray
            scaler momentum transfer Q
        dq_array: ndarray
            scaler momentum transfer (Q) resolution
        iq_array: ndarray
            I(Q)
        sigmaq_array: ndarray
            sigma I(Q)
        bin_centers: numpy.ndarray
            bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
        bin_edges: numpy.ndarray
            bin edges
        Returns
        -------
        IofQ
            named tuple for Q, dQ, binned I(Q), binned sigma_I(Q)
        """
        # check input
        assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

        # Flatten input data to 1D
        q_array = IofQCalculator.flatten(q_array)
        dq_array = IofQCalculator.flatten(dq_array)
        iq_array = IofQCalculator.flatten(iq_array)
        sigmaq_array = IofQCalculator.flatten(sigmaq_array)

        # calculate 1/sigma^2 for multiple uses
        invert_sigma2_array = 1./(sigmaq_array**2)

        print('I(Q) array:\n', iq_array)
        print('invert_sigma2_array\n', invert_sigma2_array)
        print('Raw raw I:\n', iq_array*invert_sigma2_array)
        print(bin_edges)

        # Counts per bin: I_{k, raw} = \sum \frac{I(i, j)}{(\sigma I(i, j))^2}
        i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array*invert_sigma2_array)

        print('[DEBUG 1] Raw: {}'.format(i_raw_array))
        for i in range(i_raw_array.shape[0]):
            print('{}    {:.7f}    {:.7f}    {:.7f}   {}'
                  ''.format(i, bin_x[i], bin_centers[i], bin_x[i+1], i_raw_array[i]))

        # Weight per bin: w_k = \sum \frac{1}{\sqrt{I(i, j)^2}
        w_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=invert_sigma2_array)

        # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{w_k}
        #       sigma = 1/sqrt(w_k)
        i_final_array = i_raw_array / w_array
        sigma_final_array = 1/np.sqrt(w_array)

        # Calculate Q resolution of binned
        # FIXME - waiting for Lisa's equations for binned q resolution
        # FIXME - this is an incorrect solution temporarily for workflow
        binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
        bin_q_resolution = binned_dq / i_raw_array

        # Get the final result
        binned_iq = IofQ(bin_centers, bin_q_resolution, i_final_array, sigma_final_array)

        return binned_iq

    @staticmethod
    def no_weight_binning(q_array, dq_array, iq_array, sigmaq_array, bin_centers, bin_edges):
        """ Bin I(Q) by given bin edges and do no-weight binning
        This method implements equation 11.34, 11.35 and 11.36 in master document.
        Parameters
        ----------
        q_array: ndarray
            scaler momentum transfer Q
        dq_array: ndarray
            scaler momentum transfer (Q) resolution
        iq_array: ndarray
            I(Q)
        sigmaq_array: ndarray
            sigma I(Q)
        bin_centers: numpy.ndarray
            bin centers. Note not all the bin center is center of bin_edge(i) and bin_edge(i+1)
        bin_edges: numpy.ndarray
            bin edges
        Returns
        -------
        IofQ
            named tuple for Q, dQ, binned I(Q), binned sigma_I(Q)
        """
        # check input
        assert bin_centers.shape[0] + 1 == bin_edges.shape[0]

        # Flatten input data to 1D
        q_array = IofQCalculator.flatten(q_array)
        dq_array = IofQCalculator.flatten(dq_array)
        iq_array = IofQCalculator.flatten(iq_array)
        sigmaq_array = IofQCalculator.flatten(sigmaq_array)

        # # calculate 1/sigma^2 for multiple uses
        # invert_sigma2_array = 1./(sigmaq_array**2)
        #
        # # DEBUG OUTPUT SESSION
        # print('I(Q) array:\n', iq_array)
        # print('invert_sigma2_array\n', invert_sigma2_array)
        # print('Raw raw I:\n', iq_array*invert_sigma2_array)
        print('BIN EDGES: {}'.format(bin_edges))
        # -------------------------------------------------

        # Number of I(q) in each target Q bin
        num_pt_array, bin_x = np.histogram(q_array, bins=bin_edges)

        # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
        i_raw_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=iq_array)
        # Square of summed uncertainties for each bin
        sigma_sqr_array, bin_x = np.histogram(q_array, bins=bin_edges, weights=sigmaq_array**2)

        # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{Nk}
        #       sigma = 1/sqrt(w_k)
        i_final_array = i_raw_array / num_pt_array
        sigma_final_array = np.sqrt(sigma_sqr_array) / num_pt_array

        # Calculate Q resolution of binned
        # FIXME - waiting for Lisa's equations for binned q resolution
        # FIXME - this is an incorrect solution temporarily for workflow
        binned_dq, bin_x = np.histogram(q_array, bins=bin_edges, weights=dq_array)
        bin_q_resolution = binned_dq / num_pt_array

        print('[DEBUG] No-weight: Index, Bin Center, Bin Left, Bin Right, No. Pt, I(raw), I, sigmaI')
        for i in range(i_raw_array.shape[0]):
            print('{} {:.7f} {:.7f} {:.7f} {} {} {:.7f} {:.7f}'
                  ''.format(i, bin_centers[i],  bin_x[i], bin_x[i+1], num_pt_array[i], int(i_raw_array[i]),
                            i_final_array[i], sigma_final_array[i]))

        # Get the final result
        binned_iq = IofQ(bin_centers, bin_q_resolution, i_final_array, sigma_final_array)

        return binned_iq

    @staticmethod
    def flatten(any_array):
        """
        If input array is not in shape (n, ), flatten it to (n, )
        Parameters
        ----------
        any_array : numpy.ndarray
            numpy array with at most 2 dimension
        Returns
        -------
        numpy.ndarray
            flattened input array
        """
        if len(any_array.shape) == 1:
            # no need to flatten
            return any_array

        return any_array.flatten()

    @staticmethod
    def determine_linear_bin_edges(q_min, q_max, bins):
        delta_q = (q_max - q_min) / bins
        bin_edges = np.arange(bins + 1).astype('float') * delta_q + q_min
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5

        return bin_centers, bin_edges

    @staticmethod
    def determine_log_bin_edges(q_min, q_max, step_per_decade):
        """

        Parameters
        ----------
        q_min
        q_max
        step_per_decade: float
            step per decade (ex. 0.1 to 1.0 is one decade); denoted as 'j' in document
        Returns
        -------

        """
        # Calculate step
        delta = np.power(10., 1. / step_per_decade)

        # Align q_min to power of 10 as q0
        q0 = np.power(10, np.floor(np.log10(q_min)))
        # print('[DEBUG OUTPUT: q0 = {}'.format(q0))

        # Determine number of bins
        num_bins = 1 + int(np.ceil(step_per_decade * np.log(q_max/q0)/np.log(10)))
        # print('[DEBUG OUTPUT: number of bins = {}'.format(num_bins))

        # Calculate bin centers
        bin_centers = np.arange(num_bins).astype('float')
        bin_centers = q0 * np.power(delta, bin_centers)
        # for i in range(num_bins):
        #     print('[DEBUG OUTPUT] Q({}) = {}'.format(i+1, bin_centers[i]))

        # Calculate bin boundaries
        delta_q_array = 2. * (delta - 1) / (delta + 1) * bin_centers
        bin_edges = np.zeros((num_bins+1,), dtype='float')
        bin_edges[1:] = bin_centers[:] + 0.5 * delta_q_array[:]
        bin_edges[0] = bin_centers[0] - 0.5 * delta_q_array[0]

        # # Big debug
        # print('[DEBUG OUTPUT] Edge from {}'.format(bin_edges[0]))
        # for i in range(99):
        #     from_left = bin_centers[i] + 0.5 * delta_q_array[i]
        #     from_right = bin_centers[i + 1] - 0.5 * delta_q_array[i + 1]
        #     diff = from_right - from_left
        #     diff2 = from_right - bin_edges[i+1]
        #     print('[DEBUG OUTPUT] From left = {}, From right = {}, Difference = {}, Calculated = {}  Diff = {}'
        #           ''.format(from_left, from_right, diff, bin_edges[i+1], diff2))
        # # END-FOR

        return bin_centers, bin_edges
