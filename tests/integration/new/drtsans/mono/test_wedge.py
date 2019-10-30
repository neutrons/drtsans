import pytest
import numpy as np
from collections import namedtuple
from mantid.simpleapi import CreateWorkspace, PlotPeakByLogValue, SaveNexus
from drtsans.settings import unique_workspace_dundername

def test_wedge():
    intensities = np.array([[1.5, 1.9, 2.2, 1.8, 1.9, 2.8, 2.2, 2.9, 2.2, 2.2, 1.5],
                            [2, 2.2, 1.8, 2, 3, 3.3, 2.3, 2.8, 1.8, 2.8, 3.2],
                            [2.5, 3.2, 2.2, 2.2, 2.5, 3.4, 3.4, 2.9, 2.2, 2.2, 3],
                            [4, 5.5, 2.3, 2.8, 3.8, 2.4, 2, 2.2, 7, 6, 3.5],
                            [5, 6, 12, 15, 4, 2.5, 2.4, 15, 10, 7, 4],
                            [6, 7, 15, 22, 31, 2, 30, 23, 13, 8, 5],
                            [5, 6, 11, 15, 2, 1.8, 2.2, 15, 10, 7, 4],
                            [4, 5.5, 4.2, 2, 2.2, 3.2, 3.2, 2.4, 7, 6, 3.5],
                            [2, 2, 3, 2.1, 4.3, 3.8, 2.8, 2.9, 4, 3.8, 3.3],
                            [1.8, 1.5, 2.2, 2.3, 2.3, 2.9, 1.9, 3, 1.8, 2.8, 3.1],
                            [2, 1.9, 1.5, 1.9, 2.2, 2.1, 2.3, 2.5, 2.1, 1.9, 1.8]],
                           dtype=np.float)

    errors = np.array([[0.7, 0.4, 0.6, 0.5, 0.5, 0.7, 0.6, 0.7, 0.5, 0.4, 0.75],
                       [0.5, 0.6, 0.6, 0.7, 0.6, 0.8, 0.5, 0.6, 0.5, 0.3, 0.25],
                       [0.4, 0.4, 0.8, 0.7, 0.8, 1, 0.9, 0.7, 0.6, 0.4, 0.25],
                       [0.2, 0.25, 0.7, 0.7, 0.9, 0.8, 1.1, 0.9, 0.5, 0.3, 0.3],
                       [0.2, 0.3, 0.4, 0.5, 0.2, 0.1, 0.2, 0.7, 0.3, 0.2, 0.4],
                       [0.15, 0.2, 0.5, 0.7, 1, 2, 1, 0.7, 0.4, 0.2, 0.3],
                       [0.2, 0.3, 0.6, 0.6, 0.5, 0.4, 0.8, 0.5, 0.4, 0.2, 0.25],
                       [0.25, 0.4, 0.5, 0.8, 1, 0.5, 1, 0.8, 0.6, 0.4, 0.35],
                       [0.7, 0.8, 0.7, 0.6, 1.1, 0.9, 0.7, 0.5, 0.4, 0.35, 0.3],
                       [0.5, 0.4, 0.5, 0.7, 0.8, 0.9, 0.6, 0.6, 0.6, 0.4, 0.3],
                       [0.5, 0.5, 0.7, 0.6, 0.7, 0.5, 0.6, 0.5, 0.4, 0.45, 0.5]],
                      dtype=np.float)

    # Create a Qx-Qy bi-dimensional grid
    TwoDGridXY = namedtuple('TwoDGridXY', 'x y')
    qy, qx = np.array(range(10, -1, -1)), np.array(range(11))
    qx0, qy0 = [np.mean(q_i) for q_i in (qx, qy)]  # average Q values on each dimension
    assert (qx0, qy0) == pytest.approx((5.0, 5.0))
    qx_qy_grid = TwoDGridXY(*np.meshgrid(qx, qy))  # bidimensional grid implemented as a namedtuple

    # Calculate Q and Phi for each pixel
    q = np.hypot(qx_qy_grid.x - qx0, qx_qy_grid.y - qy0)
    phi = np.degrees(np.arctan2(qx_qy_grid.y - qy0, qx_qy_grid.x - qx0))
    negative_phi = phi < 0
    phi += 360 * negative_phi.astype(float)

    # Create 2D histogram in Q and Phi using the (qx, qy, intensity) data from all the pixels.
    q_bins, phi_n_bins = np.arange(1, 9), 80
    intensities_q_phi, _, phi_bins = np.histogram2d(q.flatten(), phi.flatten(), bins=[q_bins, phi_n_bins],
                                                    weights=intensities.flatten())
    # Also create 2D histogram for errors, which are added in quadrature
    errors_q_phi, _, _ = np.histogram2d(q.flatten(), phi.flatten(), bins=[q_bins, phi_n_bins],
                                        weights=errors.flatten()**2)
    errors_q_phi = np.sqrt(errors_q_phi)

    # Create a Q-Phi bi-dimensional grid using the center points for the Q and Phi bins
    q_centers, phi_centers = (q_bins[0:-1] + q_bins[1:]) / 2, (phi_bins[0:-1] + phi_bins[1:]) / 2
    TwoDGridQPhi = namedtuple('TwoDGridQPhi', 'q phi')
    q_phi_grid = TwoDGridQPhi(*np.meshgrid(q_centers, phi_centers))

    # Fitting each Q-slice
    q_phi_workspace = unique_workspace_dundername()
    CreateWorkspace(DataX=phi_centers, DataY=intensities_q_phi, DataE=errors_q_phi, Nspec=len(q_centers),
                    UnitX='Degrees', OutputWorkspace=q_phi_workspace,
                    VerticalAxisUnit='MomentumTransfer', VerticalAxisValues=[str(q) for q in q_centers])
    SaveNexus(InputWorkspace=q_phi_workspace, Filename='/tmp/junk.nxs')
    function = 'name=FlatBackground,A0=0;name=Gaussian,Height=0,PeakCentre=180,Sigma=10'

if __name__ == '__main__':
    pytest.main([__file__])
