from __future__ import print_function

import numpy as np
from scipy import stats

from mantid.simpleapi import CreateEmptyTableWorkspace, CreateWorkspace
from ornl.sans.detector import Component
from ornl.sans.hfir import resolution

# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


def log_space(start, stop, num=100):
    """
    Return numbers spaced evenly on a log scale.

    In linear space, the sequence starts at base ** start (base to the power
    of start) and ends with base ** stop.

    > log_space(0.001, 0.6)
      [0.001,
        0.0010667487284100271,
        0.0011379528495644108,
        0.0012139097552634022,
        0.0012949366878317627,
        0.001381372065116025,
        (..)
        0.40717250825126383,
        0.43435075542055696,
        0.4633431160288143,
        0.49427067984127726,
        0.5272626192110421,
        0.5624567285815201,
        0.6]

    Parameters
    ----------
    start : float
        The starting value of the sequence.

    stop : float
        The final value of the sequence.
        Num values are spaced over the interval in log-space, of which all
        (a sequence of length num) are returned.

    num : integer, optional
        Number of samples to generate. Default is 100.

    Returns
    -------
    list
        An array with the space
    """

    return list(np.geomspace(start, stop, num=num, endpoint=True))


class MomentumTransfer:

    '''
    Momentum Transfer class
    Olds arrays for qx, qy, dqx, dqy, i, i_sigma
    q1d and q2d are calculated from there arrays
    '''

    qx, qy, dqx, dqy, i, i_sigma = np.empty(0), np.empty(0), np.empty(0), \
        np.empty(0), np.empty(0), np.empty(0)
    prefix = None
    component = None

    def __init__(self, input_workspace=None, component_name="detector1",
                 out_ws_prefix="ws"):
        self.prefix = out_ws_prefix
        if input_workspace is not None:
            self.component = Component(input_workspace, component_name)
            self._initialize_qs(input_workspace)

    def _remove_monitors(self, data):
        return data[self.component.first_index:
                    self.component.first_index + self.component.dims]

    def _mask_pixels(self, data, masked_pixels):
        return np.ma.MaskedArray(data, masked_pixels, dtype=np.float,
                                 fill_value=np.nan)

    def _initialize_qs(self, input_workspace):

        assert input_workspace.blocksize() == 1  # sanity check: only 1 bin

        masked_pixels = self.component.masked_ws_indices()

        # 1D arrays
        self.qx, self.qy, self.dqx, self.dqy = \
            resolution.q_resolution_per_pixel(input_workspace)
        self.i = input_workspace.extractY().ravel()
        self.i_sigma = input_workspace.extractE().ravel()

        # Get rid of the monitors; from 49154 to 49152 spectra
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._remove_monitors(d) for d in [
                self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]]

        # Create numpy mask arrays with the masked pixels
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._mask_pixels(d, masked_pixels) for d in [
                self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]]

    def _create_table_ws(self):
        table_iq = CreateEmptyTableWorkspace(
            OutputWorkspace=self.prefix+"_iqxqy_table")
        table_iq.addColumn(type="float", name="Qx")
        table_iq.addColumn(type="float", name="Qy")
        table_iq.addColumn(type="float", name="dQx")
        table_iq.addColumn(type="float", name="dQy")
        table_iq.addColumn(type="float", name="I")
        table_iq.addColumn(type="float", name="Sigma(I)")

        for qx_i, qy_i, dqx_i, dqy_i, i_i, i_sigma_i in zip(
                self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma):
            nextRow = {
                'Qx': qx_i,
                'Qy': qy_i,
                'dQx': dqx_i,
                'dQy': dqy_i,
                "I": i_i,
                "Sigma(I)": i_sigma_i,
            }
            table_iq.addRow(nextRow)
        return table_iq

    def q2d(self):
        return self._create_table_ws()

    def bin_into_q2d(self, bins=None):
        """Bin the data into Q 2D for visualization only! No dqx/y data!!!

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

        """

        if bins is None:
            bins = [self.component.dim_x, self.component.dim_y]

        # Number of bins in Qx Qy is the number of pixels in X and Y
        counts_qx_qy_weights, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins)
        counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins, weights=self.i)
        counts_qx_qy /= counts_qx_qy_weights

        # TODO: Error propagation is probably wrong
        counts_dqx_dqy_weights, _, _ = np.histogram2d(
            self.dqx, self.dqy, bins=bins)
        counts_dqx_dqy, _, _ = np.histogram2d(
            self.dqx, self.dqy, bins=bins, weights=self.i_sigma)
        counts_dqx_dqy /= counts_dqx_dqy_weights

        qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0

        # When doing histogram2d the masks are gone
        iqxqy_ws = CreateWorkspace(
            DataX=np.tile(qx_bin_edges, len(qy_bin_centers)),
            DataY=counts_qx_qy.T,
            DataE=counts_dqx_dqy.T,
            NSpec=len(qy_bin_centers),
            UnitX='MomentumTransfer',
            VerticalAxisUnit='MomentumTransfer',
            VerticalAxisValues=qy_bin_centers,
            OutputWorkspace=self.prefix+"_iqxqy"
        )

        return iqxqy_ws.name(), iqxqy_ws

    @staticmethod
    def _bin_into_q1d(q, dq, i, i_sigma, prefix, bins=100, statistic='mean'):

        intensity_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, i, statistic=statistic, bins=bins)

        sigma_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, i_sigma, statistic=lambda array_1d: np.sqrt(
                np.sum(np.square(array_1d))) / len(array_1d), bins=bins)

        _, dq_bin_edges, _ = \
            stats.binned_statistic(dq, i, statistic=statistic, bins=bins)

        dq_bin_centers = (dq_bin_edges[1:] + dq_bin_edges[:-1]) / 2.0

        iq = CreateWorkspace(
            DataX=np.array([q_bin_edges]),
            DataY=np.array([intensity_statistic]),
            DataE=np.array([sigma_statistic]),
            Dx=dq_bin_centers,  # bin centers!!
            NSpec=1,
            UnitX='MomentumTransfer',
            YUnitLabel='Counts',
            OutputWorkspace=prefix+"_iq"
        )

        return iq.name(), iq

    def bin_into_q1d(self, bins=100, statistic='mean'):
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

        return MomentumTransfer._bin_into_q1d(q, dq, self.i, self.i_sigma,
                                              self.prefix, bins, statistic)

    def bin_wedge_into_q1d(self, phi_0=0, phi_aperture=30,
                           bins=100, statistic='mean'):
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

        # Let's work in radians
        phi_0_rad = np.deg2rad(phi_0)
        phi_aperture_rad = np.deg2rad(phi_aperture)

        phi_aperture_min_rad = phi_0_rad - phi_aperture_rad/2
        phi_aperture_max_rad = phi_0_rad + phi_aperture_rad/2
        # opposite 180 degrees apart
        phi_aperture_min_pi_rad = phi_aperture_min_rad + np.pi
        phi_aperture_max_pi_rad = phi_aperture_max_rad + np.pi

        condition1 = (angle_rad > phi_aperture_min_rad) & \
            (angle_rad < phi_aperture_max_rad)
        # make angle > np.pi varying between np.pi and 2*np.pi, rather than the
        # initial -np.pi to np.pi
        angle_rad[angle_rad < 0] = 2*np.pi + angle_rad[angle_rad < 0]
        condition2 = (angle_rad > phi_aperture_min_pi_rad) & \
            (angle_rad < phi_aperture_max_pi_rad)
        # True where the wedge is located, otherwise false
        condition = condition1 | condition2

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        q = q[condition]
        i = self.i[condition]
        i_sigma = self.i_sigma[condition]

        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))
        dq = dq[condition]

        return MomentumTransfer._bin_into_q1d(
            q, dq, i, i_sigma, self.prefix+"_wedge_", bins, statistic)
