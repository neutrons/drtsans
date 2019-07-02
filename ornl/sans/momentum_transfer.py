from __future__ import print_function

import numpy as np
from scipy import stats

from mantid.simpleapi import CreateWorkspace, GroupWorkspaces
from ornl.sans.detector import Component
from ornl.sans.hfir import resolution


def bin_into_q2d(ws, component_name="detector1", out_ws_prefix="ws"):
    """Bin the data into Q 2D

    TODO: This needs refactoring. Too long too complicated.

    Parameters
    ----------
    ws : MatrixWorkspace
        This must be single bin.
        In the case of EQSANS with multiple bins this must be called in a
        `for` loop. Using multiprocessing.map` is another solution.

    component_name : str, optional
        The detector name. Either 'detector1' or `wing_detector`,
        by default "detector1"
    out_ws_prefix : str, optional
        The output workspace is prefix by this, by default "ws".

    Returns
    -------
    WorkSpaceGroup
        A WorkSpaceGroup containing tuples (ws name, ws) of 3 workspaces.
        To date Mantid has no option to add uncertainties to the axis.
        The resolution dx,dy is set in separated workspaces suffixed by
        `_dqx` and `_dqy`.

    """
    # flake8: noqa E712

    assert ws.blocksize() == 1  # sanity check: only 1 bin

    det = Component(ws, component_name)
    masked_pixels = det.masked_ws_indices()

    # Note that the detctor is read from the lower left corner, then up
    # the tube is read first from the bottom (256)
    i = ws.extractY()  # i.shape == (49154, 1)
    i_sigma = ws.extractE()

    # All returned arrays have the same shape == (49154,)
    qx, qy, dqx, dqy = resolution.q_resolution_per_pixel(ws)

    # Get rid of the monitors; from 49154 to 49152 spectra
    qx = qx[det.first_index:det.first_index + det.dims]
    qy = qy[det.first_index:det.first_index + det.dims]
    dqx = dqx[det.first_index:det.first_index + det.dims]
    dqy = dqy[det.first_index:det.first_index + det.dims]
    i = i[det.first_index:det.first_index + det.dims]
    i_sigma = i_sigma[det.first_index:det.first_index + det.dims]
    # get rid of the original bins: from shape == (49152, 1) to (49152,)
    i, i_sigma = i[:, 0], i_sigma[:, 0]

    # Create numpy mask arrays with the masked pixels
    qx = np.ma.MaskedArray(qx, masked_pixels, dtype=np.float,
                           fill_value=np.nan)
    qy = np.ma.MaskedArray(qy, masked_pixels, dtype=np.float,
                           fill_value=np.nan)
    dqx = np.ma.MaskedArray(
        dqx, masked_pixels, dtype=np.float, fill_value=np.nan)
    dqy = np.ma.MaskedArray(
        dqy, masked_pixels, dtype=np.float, fill_value=np.nan)
    i = np.ma.MaskedArray(i, masked_pixels, dtype=np.float, fill_value=np.nan)
    i_sigma = np.ma.MaskedArray(
        i_sigma, masked_pixels, dtype=np.float, fill_value=np.nan)

    # Number of bins in Qx Qy is the number of pixels in X and Y
    counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(
        qx, qy, bins=[det.dim_x, det.dim_y], weights=i
    )
    counts_dqx_dqy, dqx_bin_edges, dqy_bin_edges = np.histogram2d(
        dqx, dqy, bins=[det.dim_x, det.dim_y], weights=i_sigma
    )
    qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0
    # qy_bin_centers.shape == dqy_bin_centers.shape == (256,)
    dqy_bin_centers = (dqy_bin_edges[1:] + dqy_bin_edges[:-1]) / 2.0

    # Grids for I, dqx, dqy
    qx_bin_edges_grid, qy_bin_centers_grid = np.meshgrid(
        qx_bin_edges, qy_bin_centers)
    # Since the tubes are read first: reshape(n-rows,n-cols)
    # We need to fill in the columns first: 'F'
    i_grid = i.reshape(det.dim_y, det.dim_x, order='F')
    i_sigma_grid = i_sigma.reshape(det.dim_y, det.dim_x, order='F')
    # np.tile(array, (n-times, 1)) - Stack vertically n-times the array
    dqx_bin_centers_grid = np.tile(dqx_bin_edges, (len(dqy_bin_centers), 1))
    # array[np.newaxis].T transforms a vector (row) into a column array
    # The tile repeat that column len(dqx_bin_edges) times
    dqy_bin_centers_grid = np.tile(
        dqy_bin_centers[np.newaxis].T, (1, len(dqx_bin_edges)))

    # Q WS
    iqxqy_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=i_grid.filled(),  # 2D: mask as np.nan
        DataE=i_sigma_grid.filled(),  # 2D: mask as np.nan
        NSpec=len(qy_bin_centers),
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_iqxqy"
    )

    dqx_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=dqx_bin_centers_grid,  # 2D
        DataE=None,  # 2D
        NSpec=len(qy_bin_centers),
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_dqx",
    )

    dqy_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=dqy_bin_centers_grid,  # 2D
        DataE=None,  # 2D
        NSpec=len(qy_bin_centers),
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_dqy",
    )

    qxqy_wss_grouped = GroupWorkspaces(
        InputWorkspaces=[iqxqy_ws, dqx_ws, dqy_ws],
        OutputWorkspace=out_ws_prefix+"_qxqy")

    return [(ws.name(), ws) for ws in qxqy_wss_grouped]


def bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, bins=100, statistic='mean',
                 out_ws_prefix="ws"):
    """ Calculates: I(Q) and Dq
    The ws_* input parameters are the output workspaces from bin_into_q2d

    TODO:
    This needs refactoring. Too long too complicated.

    Parameters
    ----------
    ws_iqxqy : Workspace2D
        I(qx, qy)
    ws_dqx : Workspace2D
        dqx(qx, qy)
    ws_dqy : Workspace2D
        dqy(qx, qy)
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

    q_bin_centers_grid = np.sqrt(
        np.square(qx_bin_centers_grid) + np.square(qy_bin_centers_t_grid))

    # Calculate I(Q) and error(I(Q))
    i = ws_iqxqy.extractY()
    sigma_i = ws_iqxqy.extractE()
    assert(q_bin_centers_grid.shape == i.shape == sigma_i.shape)

    intensity_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
        q_bin_centers_grid.ravel(), i.ravel(), statistic=statistic, bins=bins)

    sigma_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
        q_bin_centers_grid.ravel(), sigma_i.ravel(),
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

    dq_intensity_statistic, dq_bin_edges, dq_binnumber = \
        stats.binned_statistic(dq_bin_centers_grid_all.ravel(), i.ravel(),
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

    intensity_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
        q_bin_centers_grid, i, statistic=statistic, bins=bins)

    sigma_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
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

    dq_intensity_statistic, dq_bin_edges, dq_binnumber = \
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
