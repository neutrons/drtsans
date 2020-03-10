'''The test data for all of these examples come from
https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/177

Much of the spreadsheet is split into smaller tests to aid in verifying the intermediate results
'''
import pytest
import numpy as np
from drtsans.dataobjects import IQazimuthal, IQmod
from drtsans import getWedgeSelection
# test internal functions
from drtsans.auto_wedge import _binInQAndAzimuthal, _fitQAndAzimuthal
from drtsans.iq import  _toQmodAndAzimuthal

def _create_2d_data():
    '''This creates the IQazimuthal data from ``Raw_Data_Anisotropic_v2_new`` the pages are

    * "Anisotropic Data - Qy vs Qx" contains the intensities with the Q values labeled in grey
    * "Anisotropic Data - EB Qy vs Qx" contains the uncertainties
    * "Anisotropic Data - Q" contains the Q-magnitude
    * "Anisotropic Data - Phi" contains the azimuthal angle
    '''
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
                       [24.5, 26.5, 38.7, 46.9, 55.7, 01.0, 54.8, 48.0, 36.1, 28.3, 22.4],
                       [22.4, 24.5, 34.6, 38.7, 20.0, 15.8, 15.5, 38.7, 31.6, 26.5, 20.0],
                       [20.0, 23.5, 15.2, 16.7, 19.5, 15.5, 14.1, 14.8, 26.5, 24.5, 18.7],
                       [15.8, 17.9, 14.8, 14.8, 15.8, 18.4, 18.4, 17.0, 14.8, 14.8, 17.3],
                       [14.1, 14.8, 13.4, 14.1, 17.3, 18.2, 15.2, 17.3, 13.4, 16.7, 17.9],
                       [12.2, 13.8, 14.8, 13.4, 13.8, 16.7, 14.8, 17.0, 14.8, 14.8, 12.2]])
    data2d = IQazimuthal(intensity=intensities, error=errors,
                         qx=np.linspace(-5., 5., 11, dtype=float),
                         qy=np.linspace(5., -5., 11, dtype=float))
    assert data2d.intensity.shape == (11 * 11,)
    return data2d


def _create_2d_histogram_data():
    '''This creates the parallel arrays of the binned data from ``Raw_Data_Anisotropic_v2_new`` the pages are

    * "Anisotropic Data - Q vs Phi" which contains the intensity of each Q column and the rows are the
      azimuthal angle values labeled with bin boundary. Both of these express bin centers.
    * "Anisotropic Data - EB Q vs Phi" which contains the uncertainties
    '''
    # numbers taken from the spreadsheet
    intensity = np.array([[3000, 2300, 1300, 800, 500, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 400, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 700, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1000, np.nan, 350, np.nan, np.nan],
                          [np.nan, 1500, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 330, np.nan],
                          [np.nan, np.nan, np.nan, 700, 380, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 310, np.nan],
                          [220, np.nan, 240, 400, np.nan, 280, 180],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 190, np.nan],
                          [np.nan, np.nan, np.nan, 290, 180, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 210, np.nan],
                          [np.nan, 320, np.nan, 300, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 280, np.nan, 250, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 190, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 230, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [180, 320, 380, 290, 210, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 230, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 430, np.nan, 190, np.nan, np.nan],
                          [np.nan, 220, np.nan, 230, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 150, np.nan],
                          [np.nan, np.nan, np.nan, 210, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 190, np.nan],
                          [200, np.nan, 200, 300, np.nan, 150, 200],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 180, np.nan],
                          [np.nan, np.nan, np.nan, 420, 200, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 200, np.nan],
                          [np.nan, 1500, np.nan, 550, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1100, np.nan, 400, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 500, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [3100, 2200, 1500, 700, 600, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 500, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1200, np.nan, 400, np.nan, np.nan],
                          [np.nan, 1500, np.nan, 550, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 250, np.nan],
                          [np.nan, np.nan, np.nan, 230, 320, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 200, np.nan],
                          [400, np.nan, 280, 220, np.nan, 220, 150],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 190, np.nan],
                          [np.nan, np.nan, np.nan, 220, 180, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 220, np.nan],
                          [np.nan, 380, np.nan, 200, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 250, np.nan, 180, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 300, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 190, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [250, 240, 340, 330, 280, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 230, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 340, np.nan, 290, np.nan, np.nan],
                          [np.nan, 200, np.nan, 300, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 220, np.nan],
                          [np.nan, np.nan, np.nan, 290, 180, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 220, np.nan],
                          [240, np.nan, 220, 220, np.nan, 280, 150],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 320, np.nan],
                          [np.nan, np.nan, np.nan, 700, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 300, np.nan],
                          [np.nan, 1500, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1000, np.nan, 350, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 700, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 400, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [3000, 2300, 1300, 800, 500, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 400, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 700, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1000, np.nan, 350, np.nan, np.nan],
                          [np.nan, 1500, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 330, np.nan],
                          [np.nan, np.nan, np.nan, 700, 380, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 310, np.nan],
                          [220, np.nan, 240, 400, np.nan, 280, 180],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 190, np.nan],
                          [np.nan, np.nan, np.nan, 290, 180, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 210, np.nan],
                          [np.nan, 320, np.nan, 300, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 280, np.nan, 250, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 190, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 230, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [180, 320, 380, 290, 210, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 230, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 430, np.nan, 190, np.nan, np.nan],
                          [np.nan, 220, np.nan, 230, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 150, np.nan],
                          [np.nan, np.nan, np.nan, 210, 220, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 190, np.nan],
                          [200, np.nan, 200, 300, np.nan, 150, 200],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 180, np.nan],
                          [np.nan, np.nan, np.nan, 420, 200, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, 200, np.nan],
                          [np.nan, 1500, np.nan, 550, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, 1100, np.nan, 400, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, 600, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, 500, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                          [3100, 2200, 1500, 700, 600, np.nan, np.nan]], dtype=float)

    error = np.array([[54.8, 48.0, 36.1, 28.3, 22.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 31.6, np.nan, 18.7, np.nan, np.nan],
                      [np.nan, 38.7, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 18.2, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, 19.5, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 17.6, np.nan],
                      [14.8, np.nan, 15.5, 20.0, np.nan, 16.7, 13.4],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.8, np.nan],
                      [np.nan, np.nan, np.nan, 17.0, 13.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.5, np.nan],
                      [np.nan, 17.9, np.nan, 17.3, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 16.7, np.nan, 15.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 13.8, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 15.2, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [13.4, 17.9, 19.5, 17.0, 14.5, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 15.2, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 20.7, np.nan, 13.8, np.nan, np.nan],
                      [np.nan, 14.8, np.nan, 15.2, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 12.2, np.nan],
                      [np.nan, np.nan, np.nan, 14.5, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.8, np.nan],
                      [14.1, np.nan, 14.1, 17.3, np.nan, 12.2, 14.1],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.4, np.nan],
                      [np.nan, np.nan, np.nan, 20.5, 14.1, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.1, np.nan],
                      [np.nan, 38.7, np.nan, 23.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 33.2, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 22.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [55.7, 46.9, 38.7, 26.5, 24.5, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 22.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 34.6, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, 38.7, np.nan, 23.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 15.8, np.nan],
                      [np.nan, np.nan, np.nan, 15.2, 17.9, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.1, np.nan],
                      [20.0, np.nan, 16.7, 14.8, np.nan, 14.8, 12.2],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.8, np.nan],
                      [np.nan, np.nan, np.nan, 14.8, 13.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.8, np.nan],
                      [np.nan, 19.5, np.nan, 14.1, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 15.8, np.nan, 13.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 17.3, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 13.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [15.8, 15.5, 18.4, 18.2, 16.7, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 15.2, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 18.4, np.nan, 17.0, np.nan, np.nan],
                      [np.nan, 14.1, np.nan, 17.3, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.8, np.nan],
                      [np.nan, np.nan, np.nan, 17.0, 13.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.8, np.nan],
                      [15.5, np.nan, 14.8, 14.8, np.nan, 16.7, 12.2],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 17.9, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 17.3, np.nan],
                      [np.nan, 38.7, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 31.6, np.nan, 18.7, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [54.8, 48.0, 36.1, 28.3, 22.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 31.6, np.nan, 18.7, np.nan, np.nan],
                      [np.nan, 38.7, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 18.2, np.nan],
                      [np.nan, np.nan, np.nan, 26.5, 19.5, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 17.6, np.nan],
                      [14.8, np.nan, 15.5, 20.0, np.nan, 16.7, 13.4],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.8, np.nan],
                      [np.nan, np.nan, np.nan, 17.0, 13.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.5, np.nan],
                      [np.nan, 17.9, np.nan, 17.3, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 16.7, np.nan, 15.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 13.8, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 15.2, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [13.4, 17.9, 19.5, 17.0, 14.5, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 15.2, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 20.7, np.nan, 13.8, np.nan, np.nan],
                      [np.nan, 14.8, np.nan, 15.2, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 12.2, np.nan],
                      [np.nan, np.nan, np.nan, 14.5, 14.8, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.8, np.nan],
                      [14.1, np.nan, 14.1, 17.3, np.nan, 12.2, 14.1],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 13.4, np.nan],
                      [np.nan, np.nan, np.nan, 20.5, 14.1, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, 14.1, np.nan],
                      [np.nan, 38.7, np.nan, 23.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, 33.2, np.nan, 20.0, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, 24.5, np.nan, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, 22.4, np.nan, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                      [55.7, 46.9, 38.7, 26.5, 24.5, np.nan, np.nan]], dtype=float)
    assert intensity.shape == error.shape

    # put together binning parameters to reproduce spreadsheet
    q_min = .5
    q_max = 7.5
    q_delta = 1.
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)

    azimuthal_delta = 5.
    azimuthal_max = 540. + azimuthal_delta
    azimuthal_bins = np.arange(-.5 * azimuthal_delta, azimuthal_max, azimuthal_delta, dtype=float)

    assert intensity.shape == (len(azimuthal_bins) - 1, len(q_bins) - 1)

    return intensity, error, azimuthal_bins, q_bins


def test_calc_qmod_and_azimuthal():
    '''Test the conversion of data into 2d arrays of qmod and azimuthal angle. The results are
    checked against "Anisotropic Data -Q" and "Anisotropic Data - Phi"
    '''
    data2d = _create_2d_data()

    # convert to q and azimuthal
    qmod, delta_qmod, azimuthal = _toQmodAndAzimuthal(data2d)
    assert qmod.shape == data2d.intensity.shape
    assert delta_qmod is None
    assert azimuthal.shape == data2d.intensity.shape

    # numbers taken from the spreadsheet
    q_exp = np.array([[7.07, 6.40, 5.83, 5.39, 5.10, 5.00, 5.10, 5.39, 5.83, 6.40, 7.07],
                      [6.40, 5.66, 5.00, 4.47, 4.12, 4.00, 4.12, 4.47, 5.00, 5.66, 6.40],
                      [5.83, 5.00, 4.24, 3.61, 3.16, 3.00, 3.16, 3.61, 4.24, 5.00, 5.83],
                      [5.39, 4.47, 3.61, 2.83, 2.24, 2.00, 2.24, 2.83, 3.61, 4.47, 5.39],
                      [5.10, 4.12, 3.16, 2.24, 1.41, 1.00, 1.41, 2.24, 3.16, 4.12, 5.10],
                      [5.00, 4.00, 3.00, 2.00, 1.00, 0.00, 1.00, 2.00, 3.00, 4.00, 5.00],
                      [5.10, 4.12, 3.16, 2.24, 1.41, 1.00, 1.41, 2.24, 3.16, 4.12, 5.10],
                      [5.39, 4.47, 3.61, 2.83, 2.24, 2.00, 2.24, 2.83, 3.61, 4.47, 5.39],
                      [5.83, 5.00, 4.24, 3.61, 3.16, 3.00, 3.16, 3.61, 4.24, 5.00, 5.83],
                      [6.40, 5.66, 5.00, 4.47, 4.12, 4.00, 4.12, 4.47, 5.00, 5.66, 6.40],
                      [7.07, 6.40, 5.83, 5.39, 5.10, 5.00, 5.10, 5.39, 5.83, 6.40, 7.07]], dtype=float)

    azimuthal_exp = np.array([[135, 129, 121, 112, 101, 90, 79, 68, 59, 51, 45],
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

    np.testing.assert_allclose(qmod, q_exp, atol=.005)
    np.testing.assert_allclose(azimuthal, azimuthal_exp, atol=.5)


def test_bin_into_q_and_azimuthal():
    '''Test binning into Q and azimuthal matches the results from "Anisotropic Data - Q vs Phi"'''
    # get the test data
    data2d = _create_2d_data()
    intensity_exp, error_exp, _, _ = _create_2d_histogram_data()

    # parameters for azimuthal
    azimuthal_delta = 5.

    # parameters for q
    q_min = .5
    q_max = 7.5
    q_delta = 1.

    # get the histogrammed data
    intensity, error, azimuthal, q = _binInQAndAzimuthal(data2d, q_min=q_min, q_max=q_max, q_delta=q_delta,
                                                         azimuthal_delta=azimuthal_delta)

    # verify generic shape and range values
    assert azimuthal.min() == -0.5 * azimuthal_delta  # phi bins are centered on [0, 5, 10, ...]
    assert azimuthal.max() == 540. + 0.5 * azimuthal_delta
    assert q.min() == q_min
    assert q.max() == q_max  # using bin boundaries

    assert intensity.shape == intensity_exp.shape
    assert error.shape == error_exp.shape
    # verify shape is consistent with histogramming
    assert intensity.shape == (len(azimuthal) - 1, len(q) - 1)

    # validate results
    for i in range(len(intensity)):  # loop over rows to make debugging easier
        msg = 'i={} | {}deg <= azimuthal < {}deg'.format(i+3, azimuthal[i], azimuthal[i+1])
        np.testing.assert_allclose(intensity[i], intensity_exp[i], atol=.05, equal_nan=True, err_msg=msg)
        np.testing.assert_allclose(error[i], error_exp[i], atol=.05, equal_nan=True, err_msg=msg)


def test_fitting():
    '''Test that the fitting generates reasonable results for fitting the peaks'''
    intensity, error, azimuthal, q = _create_2d_histogram_data()
    # this calling forces there to be two found peaks
    center_list, fwhm_list = _fitQAndAzimuthal(intensity, error, azimuthal, q,
                                               signal_to_noise_min=2.0,
                                               azimuthal_start=110.,
                                               maxchisq=1000.)

    assert center_list[0] == pytest.approx(180., abs=1.)
    assert center_list[1] == pytest.approx(360., abs=1.)
    assert fwhm_list[0] == pytest.approx(fwhm_list[1], abs=2.)


def test_integration():
    '''Test the full workflow of the algorithm'''
    data2d = _create_2d_data()
    # parameters for azimuthal
    azimuthal_delta = 5.

    # parameters for q
    q_min = .5
    q_max = 7.5
    q_delta = 1.

    # run the function
    mins_and_maxes = getWedgeSelection(data2d, q_min=q_min, q_max=q_max, q_delta=q_delta,
                                       azimuthal_delta=azimuthal_delta)

    # tests
    assert len(mins_and_maxes) == 4
    for min_val, max_val in mins_and_maxes:
        print(min_val, '<', max_val)
        assert -90. < min_val < 270.
        assert -90. < max_val < 270.

    # first peak
    assert 0.5 * (mins_and_maxes[0][0] + mins_and_maxes[0][1]) == pytest.approx(180., abs=1.), \
        'First center is at 180.'
    assert mins_and_maxes[0][0] == pytest.approx(169., abs=.5)
    assert mins_and_maxes[0][1] == pytest.approx(192., abs=.5)

    # first background - the extra 360 is to get around the circle
    assert 0.5 * (mins_and_maxes[1][0] + mins_and_maxes[1][1] + 360) == pytest.approx(270., abs=1.2), \
        'Second center is at 270.'
    assert mins_and_maxes[1][0] == pytest.approx(249., abs=.5)
    assert mins_and_maxes[1][1] == pytest.approx(-70., abs=.5)

    # second peak
    assert 0.5 * (mins_and_maxes[2][0] + mins_and_maxes[2][1]) == pytest.approx(0., abs=1.), 'Third center is at 0.'
    assert mins_and_maxes[2][0] == pytest.approx(-11., abs=.5)
    assert mins_and_maxes[2][1] == pytest.approx(12., abs=.5)

    # second background
    assert 0.5 * (mins_and_maxes[3][0] + mins_and_maxes[3][1]) == pytest.approx(90., abs=2.), 'Forth center is at 90.'
    assert mins_and_maxes[3][0] == pytest.approx(71., abs=.5)
    assert mins_and_maxes[3][1] == pytest.approx(112., abs=.5)


if __name__ == '__main__':
    pytest.main([__file__])
