from string import Template

import numpy as np
from scipy import stats

import mantid
from mantid.simpleapi import CreateEmptyTableWorkspace, CreateWorkspace
from mantid.api import AnalysisDataService
from drtsans.momentum_transfer_factory import calculate_q_dq
from drtsans.detector import Component

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


def bin_into_q1d(wl_ws, q_bin, statistic):
    """Binning a workspace (in wave length unit) to Q (1D scaler)
    :param wl_ws: List of workspaces (names) in binned wave length space
    :param q_bin: range and bi size for Q
    :param statistic:
    :return:
    """
    calculator = IofQCalculator(wl_ws)

    return calculator.bin_into_q1d(q_bin, statistic)


def bin_into_q2d(wl_ws, bins, suffix):
    """ Bin the data into Q in 2D (Qx, Qy)
    :param wl_ws: List of workspaces (names) in binned wave length space
    :param bins: Iterable for range and bin size of Qx and Qy
    :param suffix: suffix for output workspace
    :return:
    """
    calculator = IofQCalculator(wl_ws)

    return calculator.bin_into_q2d(bins=bins, suffix=suffix)


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

    def __init__(self, input_workspace=None, component_name="detector1", out_ws_prefix="ws"):
        """

        Parameters
        ----------
        input_workspace : String or mantid.dataobjects.Workspace2D
            Name or reference of Workspace in unit of wave length
        component_name
        out_ws_prefix : String
            Prefix for output workspace such as ....
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
        self._q_dq = calculate_q_dq(input_workspace)

        # Extract I(Q) from input workspace: N x M array
        # such that self._i_q[n, m] corresponds to self._q_dq[n, m]
        self._i_q = input_workspace.extractY().ravel()
        self._i_q_sigma = input_workspace.extractE().ravel()

        # Mask monitors and detectors pixels
        masked_pixels = component.masked_ws_indices()
        monitor_pixels = component.monitor_indices()  # TODO FIXME - check monitor_indices implemented?
        self._mask_pixels(masked_pixels, monitor_pixels)

        return

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
        self._i_q = np.ma.MaskedArray(self._i_q,
                                      masked_pixels,
                                      dtype=np.float,
                                      fill_value=np.nan)

        # Mask sigma I(Q)
        self._i_q_sigma = np.ma.MaskedArray(self._i_q_sigma,
                                            masked_pixels,
                                            dtype=np.float,
                                            fill_value=np.nan)

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

    def bin_into_q1d(self, bins=100, statistic='mean', suffix="_iq"):
        """ Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------

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

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))

        return IofQCalculator._bin_intensity_into_q1d(q, dq, self.i, self.i_sigma,
                                                      self.prefix, bins, statistic,
                                                      suffix)

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
