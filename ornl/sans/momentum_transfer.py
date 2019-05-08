from __future__ import print_function

import numpy as np
from scipy import stats

from mantid.simpleapi import CreateWorkspace, GroupWorkspaces
from ornl.sans.hfir import resolution
from ornl.sans.detector import Detector

'''
LoadHFIRSANS(Filename='/home/rhf/git/sans-rewrite/data/new/ornl/sans/hfir/biosans/BioSANS_exp327_scan0014_0001.xml', OutputWorkspace='biosans')
LoadHFIRSANS(Filename='/home/rhf/git/sans-rewrite/data/new/ornl/sans/hfir/gpsans/CG2_exp206_scan0017_0001.xml', OutputWorkspace='flood')
import os
os.chdir("/home/rhf/git/sans-rewrite")
ws = mtd['flood']
from ornl.sans.hfir.momentum_transfer import bin_into_q2d, bin_into_q1d
bin_into_q2d(ws)
bin_into_q1d(mtd['ws_iqxqy'], mtd['ws_dqx'], mtd['ws_dqy'])
out_ws_prefix="ws"
import numpy as np
from ornl.sans.hfir import resolution
'''

'''

# Run in //

from multiprocessing import Pool, cpu_count
from mantid.simpleapi import CreateWorkspace, AnalysisDataService

def f(_):
    dataX = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    dataY = [1,2,3,4,5,6,7,8,9,10,11,12]
    dataE = [1,2,3,4,5,6,7,8,9,10,11,12]
    dX = [1,2,3,4,5,6,7,8,9,10,11,12]
    ws = CreateWorkspace(
        DataX=dataX, DataY=dataY, DataE=dataE, NSpec=4, UnitX="Wavelength",
        Dx=dX)
    return ws

p = Pool(cpu_count())
wss_names = ["ws_{}".format(suffix) for suffix in range(10)]
wss = p.map(f, wss_names)
[AnalysisDataService.add(name, ws) for name, ws in zip(wss_names, wss)]



'''


def bin_into_q2d(ws, component_name="detector1", out_ws_prefix="ws"):

    det = Detector(ws)
    X_SIZE_DET, Y_SIZE_DET = det.get_detector_dimensions(component_name)
    N_MONITORS = det.get_number_of_monitors()

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


def bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, bins=100, statistic='mean', out_ws_prefix="ws"):
    '''
    Calulates:
    I(Q) and Dq

    bins : int or sequence of scalars, optional
    statistic : string or callable, optional: sum, mean, median
    '''

    #
    # Calculate Q
    qx_bin_edges_grid = ws_iqxqy.extractX()
    qy_bin_centers = ws_iqxqy.getAxis(1).extractValues()
    # qy_bin_centers.shape == (256,)

    # Qx
    qx_bin_centers_grid = (qx_bin_edges_grid[:, 1:] + qx_bin_edges_grid[:, :-1]) / 2.0
    # qx_bin_centers_grid.shape == (256, 192)
    
    # Qy
    qy_bin_centers_t = np.transpose([qy_bin_centers])
    qy_bin_centers_t_grid = np.tile(qy_bin_centers_t, qx_bin_edges_grid.shape[1]-1)
    # qy_bin_centers_t_grid.shape == (256, 192)

    q_bin_centers_grid = np.sqrt(np.square(qx_bin_centers_grid) + np.square(qy_bin_centers_t_grid))

    #
    # Calculate I(Q) and error(I(Q))
    i = ws_iqxqy.extractY()
    sigma_i = ws_iqxqy.extractE()
    assert(q_bin_centers_grid.shape == i.shape == sigma_i.shape)  # sanity check

    intensity_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
        q_bin_centers_grid.ravel(), i.ravel(), statistic=statistic, bins=bins)

    sigma_statistic, q_bin_edges, q_binnumber = stats.binned_statistic(
        q_bin_centers_grid.ravel(), sigma_i.ravel(),
        statistic=lambda array_1d: np.sqrt(
            np.sum(np.square(array_1d))) / len(array_1d),
        bins=bins)
    
    #
    # Calculate dq from dqx dqy
    dqx_bin_centers_grid = ws_dqx.extractY()
    dqy_bin_centers_grid = ws_dqy.extractY()

    # Bin centres in y edges in x
    dq_bin_centers_grid = np.sqrt(np.square(dqx_bin_centers_grid) + np.square(dqy_bin_centers_grid))
    # get all to centres
    dq_bin_centers_grid_all = (dq_bin_centers_grid[:,1:] + dq_bin_centers_grid[:,:-1]) / 2.0

    dq_intensity_statistic, dq_bin_edges, dq_binnumber = stats.binned_statistic(
        dq_bin_centers_grid_all.ravel(), i.ravel(), statistic=statistic, bins=bins)

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

    return iq
