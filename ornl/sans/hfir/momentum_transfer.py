from __future__ import print_function
from ornl.sans.hfir import resolution
from mantid.simpleapi import CreateWorkspace, GroupWorkspaces
import numpy as np


def bin_into_q2d(ws):

    # TODO: get this somewhere
    X_SIZE_DET = 192
    Y_SIZE_DET = 256
    N_MONITORS = 2

    # Get WS data
    i = ws.extractY()
    i_sigma = ws.extractE()

    # All returned arrays have the same shape. They are a list of spectrum = 1D
    qx, qy, dqx, dqy = resolution.q_resolution_per_pixel(ws)
    qx, qy, dqx, dqy = qx[N_MONITORS:], qy[N_MONITORS:], dqx[N_MONITORS:], \
        dqy[N_MONITORS:]

    qx_bin_centers = np.linspace(np.abs(qx).min(), np.abs(qx).max(), X_SIZE_DET)
    qy_bin_centers = np.linspace(np.abs(qy).min(), np.abs(qy).max(), Y_SIZE_DET)

    qx_bin_dist = abs(qx_bin_centers[1]-qx_bin_centers[0])
    qx_bin_edges = qx_bin_centers - qx_bin_dist/2
    qx_bin_edges = np.append(qx_bin_edges, qx_bin_centers[-1] + qx_bin_dist/2)

    qx_histo, _ = np.histogram(qx, bins=qx_bin_edges)

    iq, _ = np.histogram(q, bins=bins)
    dqx, _ = np.histogram(q, bins=bins, weights=_dqx)
    dqy, _ = np.histogram(q, bins=bins, weights=_dqy)

    # Q WS
    iqxqy_ws2d = CreateWorkspace(
        DataX=qx_bin_grid,  # 2D
        DataY=i,  # 2D
        DataE=i_sigma,  # 2D
        NSpec=256,
        UnitX='MomentumTransfer',
        VerticalAxisUnit='MomentumTransfer',
        VerticalAxisValues=qy_bin_centers,  # 1D
    )


def bin_into_q1d(ws):
    pass
