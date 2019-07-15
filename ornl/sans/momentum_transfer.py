from __future__ import print_function

import numpy as np
from scipy import stats

from mantid.simpleapi import CreateWorkspace, GroupWorkspaces, CreateEmptyTableWorkspace
from ornl.sans.detector import Component
from ornl.sans.hfir import resolution

class MomentumTransfer:

    qx, qy, dqx, dqy, i, i_sigma = None, None, None, None, None, None
    prefix = None
    component = None

    def __init__(self, input_workspace, component_name="detector1",
                 out_ws_prefix="ws"):
        self.prexix = out_ws_prefix
        self.component = Component(input_workspace, component_name)
        self._initialize_qs(input_workspace)

    def _initialize_qs(self, input_workspace):
        
        assert input_workspace.blocksize() == 1  # sanity check: only 1 bin

        masked_pixels = self.component.masked_ws_indices()

        # 1D arrays
        self.qx, self.qy, self.dqx, self.dqy = resolution.q_resolution_per_pixel(
            input_workspace)
        self.i = input_workspace.extractY().ravel().ravel()
        self.i_sigma = input_workspace.extractE().ravel()

        # Get rid of the monitors; from 49154 to 49152 spectra
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._remove_monitors(d) for d in [
                self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]]

        # Create numpy mask arrays with the masked pixels
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._mask_pixels(d, masked_pixels) for d in [
                self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]]

    def q2d(self):
        self._create_table_ws()


    def _remove_monitors(self, data):
        return data[self.component.first_index:
            self.component.first_index + self.component.dims]


    def _mask_pixels(self, data, masked_pixels):
        return np.ma.MaskedArray(data, masked_pixels, dtype=np.float,
                                fill_value=np.nan)


    def _create_table_ws(self):
        table_iq = CreateEmptyTableWorkspace(OutputWorkspace=self.prefix+"_table")
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


    def bin_into_q2d(self, bins=None):
        """Bin the data into Q 2D

        Parameters
        ----------

        bins : int or array_like or [int, int] or [array, array], optional
            The bin specification:
                - If int, the number of bins for the two dimensions (nx=ny=bins).
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
        WorkSpaceGroup
            A WorkSpaceGroup containing tuples (ws name, ws) of 3 workspaces.
            To date Mantid has no option to add uncertainties to the axis.
            The resolution dx,dy is set in separated workspaces suffixed by
            `_dqx` and `_dqy`.

        """
        
        if bins is None:
            bins = [self.component.dim_x, self.component.dim_y]
        
        # Number of bins in Qx Qy is the number of pixels in X and Y
        counts_qx_qy_weights, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins)
        counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins, weights=self.i)
        counts_qx_qy /= counts_qx_qy_weights

        # TODO:,Error propagation is wrong
        counts_dqx_dqy_weights, _, _ = np.histogram2d(
            self.dqx, self.dqy, bins=bins)
        counts_dqx_dqy, _, _ = np.histogram2d(
            self.dqx, self.dqy, bins=bins, weights=self.i_sigma)
        counts_dqx_dqy /= counts_dqx_dqy_weights

        qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0

        # Q WS
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


    def bin_into_q1d(self, bins=100, statistic='mean'):
        """ Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d



        Parameters
        ----------

        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (10 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
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

        intensity_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, self.i, statistic=statistic, bins=bins)

        sigma_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, self.i_sigma, statistic=lambda array_1d: np.sqrt(
                np.sum(np.square(array_1d))) / len(array_1d), bins=bins)

        # Bin centres in y edges in x
        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))

        _, dq_bin_edges, _ = \
            stats.binned_statistic(dq, self.i, statistic=statistic, bins=bins)

        dq_bin_centers = (dq_bin_edges[1:] + dq_bin_edges[:-1]) / 2.0

        iq = CreateWorkspace(
            DataX=np.array([q_bin_edges]),
            DataY=np.array([intensity_statistic]),
            DataE=np.array([sigma_statistic]),
            Dx=dq_bin_centers,  # bin centers!!
            NSpec=1,
            UnitX='MomentumTransfer',
            YUnitLabel='Counts',
            OutputWorkspace=self.prefix+"_iq"
        )

        return iq.name(), iq


def bin_wedge_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, phi_0=0, phi_aperture=30,
                       bins=100, statistic='mean', out_ws_prefix="ws_wedge"):
    """
    Wedge calculation and integration

    Calculates: I(Q) and Dq
    The ws_* input parameters are the output workspaces from bin_into_q2d

    TODO:
    The code is almost copy paste bin_into_q1d :(
    Ideally one would call inside this function bin_into_q1d but mantid
    does not allowing masking without going through a for cycle and mask
    by workspace index!

    Parameters
    ----------
    ws_iqxqy : Workspace2D
        I(qx, qy)
    ws_dqx : Workspace2D
        dqx(qx, qy)
    ws_dqy : Workspace2D
        dqy(qx, qy)
    phi_0 : int, optional
        Where to start the wedge, by default 0
    phi_aperture : int, optional
        Aperture of the wedge, by default 30
    bins : int or sequence of scalars, optional
        See `scipy.stats.binned_statistic`.
        If `bins` is an int, it defines the number of equal-width bins in the
        given range (10 by default).  If `bins` is a sequence, it defines the
        bin edges, including the rightmost edge, allowing for non-uniform bin
        widths.  Values in `x` that are smaller than lowest bin edge are
        assigned to bin number 0, values beyond the highest bin are assigned to
        ``bins[-1]``.  If the bin edges are specified, the number of bins will
        be, (nx = len(bins)-1).
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
        The prefix of the workspace created in Mantid, by default "ws_wedge"

    Returns
    -------
    tuple
         tuple (workspace name, workspace)

    """

    # Calculate Q: qy_bin_centers.shape == (256,)
    qx_bin_edges_grid = ws_iqxqy.extractX()
    qy_bin_centers = ws_iqxqy.getAxis(1).extractValues()

    # Qx: qx_bin_centers_grid.shape == (256, 192)
    qx_bin_centers_grid = (
        qx_bin_edges_grid[:, 1:] + qx_bin_edges_grid[:, :-1]) / 2.0

    # Qy: qy_bin_centers_t_grid.shape == (256, 192)
    qy_bin_centers_t = qy_bin_centers[np.newaxis].T
    qy_bin_centers_t_grid = np.tile(
        qy_bin_centers_t, qx_bin_edges_grid.shape[1]-1)

    # Q
    q_bin_centers_grid = np.sqrt(
        np.square(qx_bin_centers_grid) + np.square(qy_bin_centers_t_grid))

    # Angle
    angle_grid = np.arctan2(qy_bin_centers_t_grid, qx_bin_centers_grid)

    # This is just to show the angle
    CreateWorkspace(
        DataX=qx_bin_edges_grid,
        DataY=angle_grid,
        DataE=np.sqrt(angle_grid),
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,
        OutputWorkspace=out_ws_prefix+"_angle",
    )

    # Let's work in radians
    phi_0 = np.deg2rad(phi_0)
    phi_aperture = np.deg2rad(phi_aperture)

    phi_aperture_min = phi_0 - phi_aperture/2
    phi_aperture_max = phi_0 + phi_aperture/2
    # opposite 180 degrees apart
    phi_aperture_min_pi = phi_aperture_min + np.pi
    phi_aperture_max_pi = phi_aperture_max + np.pi

    condition1 = (angle_grid > phi_aperture_min) & \
        (angle_grid < phi_aperture_max)
    # make angle > np.pi varying between np.pi and 2*np.pi, rather than the
    # initial -np.pi to np.pi
    angle_grid[angle_grid < 0] = 2*np.pi + angle_grid[angle_grid < 0]
    condition2 = (angle_grid > phi_aperture_min_pi) & \
        (angle_grid < phi_aperture_max_pi)
    # 2D Array: True where the wedge is, otherwise false
    condition = condition1 | condition2

    # This is just to show the condition
    CreateWorkspace(
        DataX=qx_bin_edges_grid,
        DataY=condition,
        DataE=np.sqrt(angle_grid),
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,
        OutputWorkspace=out_ws_prefix+"_condition",
    )

    # This transforms a 2D into 1D array
    q_bin_centers_grid = q_bin_centers_grid[condition]

    # Calculate I(Q) and error(I(Q))
    i = ws_iqxqy.extractY()
    sigma_i = ws_iqxqy.extractE()

    i = i[condition]
    sigma_i = sigma_i[condition]

    assert(q_bin_centers_grid.shape == i.shape == sigma_i.shape)

    intensity_statistic, q_bin_edges, _ = stats.binned_statistic(
        q_bin_centers_grid, i, statistic=statistic, bins=bins)

    sigma_statistic, q_bin_edges, _ = stats.binned_statistic(
        q_bin_centers_grid, sigma_i,
        statistic=lambda array_1d: np.sqrt(
            np.sum(np.square(array_1d))) / len(array_1d), bins=bins)

    # Calculate dq from dqx dqy
    dqx_bin_centers_grid = ws_dqx.extractY()
    dqy_bin_centers_grid = ws_dqy.extractY()

    # Bin centres in y edges in x
    dq_bin_centers_grid = np.sqrt(
        np.square(dqx_bin_centers_grid) + np.square(dqy_bin_centers_grid))
    # get all to centres
    dq_bin_centers_grid_all = (
        dq_bin_centers_grid[:, 1:] + dq_bin_centers_grid[:, :-1]) / 2.0

    dq_bin_centers_grid_all = dq_bin_centers_grid_all[condition]

    _, dq_bin_edges, _ = \
        stats.binned_statistic(dq_bin_centers_grid_all, i,
                               statistic=statistic, bins=bins)

    dq_bin_centers = (dq_bin_edges[1:] + dq_bin_edges[:-1]) / 2.0

    iq = CreateWorkspace(
        DataX=np.array([q_bin_edges]),
        DataY=np.array([intensity_statistic]),
        DataE=np.array([sigma_statistic]),
        Dx=dq_bin_centers,  # bin centers!!
        NSpec=1,
        UnitX='MomentumTransfer',
        YUnitLabel='Counts',
        OutputWorkspace=out_ws_prefix+"_iq"
    )
    return iq.name(), iq

