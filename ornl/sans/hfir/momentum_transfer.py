from __future__ import print_function
from ornl.sans.hfir import resolution
from mantid.simpleapi import CreateWorkspace, GroupWorkspaces
import numpy as np

'''
LoadHFIRSANS(Filename='/home/rhf/git/sans-rewrite/data/new/ornl/sans/hfir/"
"gpsans/CG2_exp206_scan0017_0001.xml', OutputWorkspace='flood')
import os
os.chdir("/home/rhf/git/sans-rewrite")
ws = mtd['flood']
from ornl.sans.hfir.momentum_transfer import bin_into_q2d
bin_into_q2d(ws)
out_ws_prefix="ws"
import numpy as np
from ornl.sans.hfir import resolution
'''


def bin_into_q2d(ws, out_ws_prefix="ws"):

    # TODO: get this somewhere
    X_SIZE_DET = 192
    Y_SIZE_DET = 256
    N_MONITORS = 2

    # Get WS data
    i = ws.extractY()
    i_sigma = ws.extractE()

    # All returned arrays have the same shape. They are a list of spectrum = 1D
    qx, qy, dqx, dqy = resolution.q_resolution_per_pixel(ws)

    # Get rid of the monitors
    qx, qy, dqx, dqy = qx[N_MONITORS:], qy[N_MONITORS:], dqx[N_MONITORS:], \
        dqy[N_MONITORS:]
    i, i_sigma = i[N_MONITORS:], i_sigma[N_MONITORS:]
    # get rid of the original bins, transform in 1D
    i, i_sigma = i[:, 0], i_sigma[:, 0]

    # Number of bins is the number of pixels in each diraction
    counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(
        qx, qy, bins=[X_SIZE_DET, Y_SIZE_DET], weights=i
    )
    counts_dqx_dqy, dqx_bin_edges, dqy_bin_edges = np.histogram2d(
        dqx, dqy, bins=[X_SIZE_DET, Y_SIZE_DET], weights=i_sigma
    )

    qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0
    dqy_bin_centers = (dqy_bin_edges[1:] + dqy_bin_edges[:-1]) / 2.0

    # Grids for I, dqx, dqy
    qx_bin_edges_grid, qy_bin_centers_grid = np.meshgrid(
        qx_bin_edges, qy_bin_centers)
    i_grid = i.reshape(X_SIZE_DET, Y_SIZE_DET)
    i_sigma_grid = i_sigma.reshape(X_SIZE_DET, Y_SIZE_DET)

    dqx_bin_centers_grid = np.tile(dqx_bin_edges, (len(dqy_bin_centers), 1))
    dqy_bin_centers_grid = np.tile(np.array([dqy_bin_centers]).transpose(),
                                   (1, len(dqx_bin_edges)))
    # Q WS
    iqxqy_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=i_grid.T,  # 2D
        DataE=i_sigma_grid.T,  # 2D
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_iqxqy"
    )

    dqx_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=dqx_bin_centers_grid,  # 2D
        DataE=None,  # 2D
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_dqx",
    )

    dqy_ws = CreateWorkspace(
        DataX=qx_bin_edges_grid,  # 2D
        DataY=dqy_bin_centers_grid,  # 2D
        DataE=None,  # 2D
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
        OutputWorkspace=out_ws_prefix+"_dqy",
    )

    qxqy_wss_grouped = GroupWorkspaces(
        InputWorkspaces=[iqxqy_ws, dqx_ws, dqy_ws],
        OutputWorkspace=out_ws_prefix+"_qxqy")

    return qxqy_wss_grouped


def bin_into_q1d(ws):
    pass
