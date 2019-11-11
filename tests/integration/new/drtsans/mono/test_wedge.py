import pytest
import numpy as np
from collections import namedtuple
from mantid.simpleapi import CreateWorkspace, SaveNexus
from drtsans.dataobjects import IQazimuthal
from drtsans.mono.wedge import _toQmodAndPhi
from drtsans.settings import unique_workspace_dundername


def _create_2d_data():
    intensities = np.array([[200, 190, 150, 190, 220, 210, 230, 250, 210, 190, 180],
                            [180, 150, 220, 230, 230, 290, 190, 300, 180, 280, 310],
                            [200, 200, 300, 210, 430, 380, 280, 290, 400, 380, 330],
                            [400, 550, 420, 200, 220, 320, 320, 240, 700, 600, 350],
                            [500, 600, 1100, 1500, 200, 180, 220, 1500, 1000, 700, 400],
                            [600, 700, 1500, 2200, 3100, 0, 3000, 2300, 1300, 800, 500],
                            [500, 600, 1200, 1500, 400, 250, 240, 1500, 1000, 700, 400],
                            [400, 550, 230, 280, 380, 240, 200, 220, 700, 600, 350],
                            [250, 320, 220, 220, 250, 340, 340, 290, 220, 220, 300],
                            [200, 220, 180, 200, 300, 330, 230, 300, 180, 280, 320],
                            [150, 190, 220, 180, 190, 280, 220, 290, 220, 220, 150]])
    errors = np.array([[14.1, 13.8, 12.2, 13.8, 14.8, 14.5, 15.2, 15.8, 14.5, 13.8, 13.4],
                       [13.4, 12.2, 14.8, 15.2, 15.2, 17.0, 13.8, 17.3, 13.4, 16.7, 17.6],
                       [14.1, 14.1, 17.3, 14.5, 20.7, 19.5, 16.7, 17.0, 20.0, 19.5, 18.2],
                       [20.0, 23.5, 20.5, 14.1, 14.8, 17.9, 17.9, 15.5, 26.5, 24.5, 18.7],
                       [22.4, 24.5, 33.2, 38.7, 14.1, 13.4, 14.8, 38.7, 31.6, 26.5, 20.0],
                       [24.5, 26.5, 38.7, 46.9, 55.7, 00.0, 54.8, 48.0, 36.1, 28.3, 22.4],
                       [22.4, 24.5, 34.6, 38.7, 20.0, 15.8, 15.5, 38.7, 31.6, 26.5, 20.0],
                       [20.0, 23.5, 15.2, 16.7, 19.5, 15.5, 14.1, 14.8, 26.5, 24.5, 18.7],
                       [15.8, 17.9, 14.8, 14.8, 15.8, 18.4, 18.4, 17.0, 14.8, 14.8, 17.3],
                       [14.1, 14.8, 13.4, 14.1, 17.3, 18.2, 15.2, 17.3, 13.4, 16.7, 17.9],
                       [12.2, 13.8, 14.8, 13.4, 13.8, 16.7, 14.8, 17.0, 14.8, 14.8, 12.2]])
    data2d = IQazimuthal(intensity=intensities, error=errors,
                         qx=np.linspace(-5., 5., 11, dtype=float),
                         qy=np.linspace(5., -5., 11, dtype=float))
    assert data2d.intensity.shape == (11, 11)
    return data2d


def test_qmod_and_phi():
    data2d = _create_2d_data()

    # convert to q and phi
    qmod, phi = _toQmodAndPhi(data2d)
    assert qmod.shape == data2d.intensity.shape
    assert phi.shape == data2d.intensity.shape

    # numbers taken from the spreadsheet
    q_exp = np.array([[7.1, 6.4, 5.8, 5.4, 5.1, 5.0, 5.1, 5.4, 5.8, 6.4, 7.1],
                      [6.4, 5.7, 5.0, 4.5, 4.1, 4.0, 4.1, 4.5, 5.0, 5.7, 6.4],
                      [5.8, 5.0, 4.2, 3.6, 3.2, 3.0, 3.2, 3.6, 4.2, 5.0, 5.8],
                      [5.4, 4.5, 3.6, 2.8, 2.2, 2.0, 2.2, 2.8, 3.6, 4.5, 5.4],
                      [5.1, 4.1, 3.2, 2.2, 1.4, 1.0, 1.4, 2.2, 3.2, 4.1, 5.1],
                      [5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                      [5.1, 4.1, 3.2, 2.2, 1.4, 1.0, 1.4, 2.2, 3.2, 4.1, 5.1],
                      [5.4, 4.5, 3.6, 2.8, 2.2, 2.0, 2.2, 2.8, 3.6, 4.5, 5.4],
                      [5.8, 5.0, 4.2, 3.6, 3.2, 3.0, 3.2, 3.6, 4.2, 5.0, 5.8],
                      [6.4, 5.7, 5.0, 4.5, 4.1, 4.0, 4.1, 4.5, 5.0, 5.7, 6.4],
                      [7.1, 6.4, 5.8, 5.4, 5.1, 5.0, 5.1, 5.4, 5.8, 6.4, 7.1]], dtype=float)
    phi_exp = np.array([[135, 129, 121, 112, 101, 90, 79, 68, 59, 51, 45],
                        [141, 135, 127, 117, 104, 90, 76, 63, 53, 45, 39],
                        [149, 143, 135, 124, 108, 90, 72, 56, 45, 37, 31],
                        [158, 153, 146, 135, 117, 90, 63, 45, 34, 27, 22],
                        [169, 166, 162, 153, 135, 90, 45, 27, 18, 14, 11],
                        [180, 180, 180, 180, 180, 0, 0, 0, 0, 0, 0],
                        [191, 194, 198, 207, 225, 270, 315, 333, 342, 346, 349],
                        [202, 207, 214, 225, 243, 270, 297, 315, 326, 333, 338],
                        [211, 217, 225, 236, 252, 270, 288, 304, 315, 323, 329],
                        [219, 225, 233, 243, 256, 270, 284, 297, 307, 315, 321],
                        [225, 231, 239, 248, 259, 270, 281, 292, 301, 309, 315]], dtype=float)
    phi = np.rad2deg(phi)  # makes comparing with spreadsheet easier

    np.testing.assert_allclose(qmod, q_exp, atol=.05)
    np.testing.assert_allclose(phi, phi_exp, atol=.5)


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
