"""The test data for all of these examples come from
https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/177

Much of the spreadsheet is split into smaller tests to aid in verifying the intermediate results
"""

import pytest
import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
from unittest.mock import MagicMock
import drtsans.auto_wedge
from matplotlib.image import imread
from drtsans.auto_wedge import (
    RING_PLOT_COLUMNS,
    RING_PLOT_PANEL_SIZE,
    AcceptanceCriteria,
    AutoWedgeError,
    RingStatus,
    SelectionFailure,
    _binInQAndAzimuthal,
    _estimatePeakParameters,
    _fitQAndAzimuthal,
)
from drtsans.dataobjects import IQazimuthal, I1DAnnular, load_iq2d_from_h5
from drtsans.determine_bins import determine_1d_linear_bins
from drtsans import getWedgeSelection
from drtsans.mono import biosans
from drtsans.plots import plot_IQazimuthal

# test internal function _toQmodAndAzimuthal as well
from drtsans.iq import (
    _toQmodAndAzimuthal,
    valid_wedge,
    BinningMethod,
    bin_intensity_into_q2d,
    select_i_of_q_by_wedge,
    bin_all,
)
from matplotlib.colors import LogNorm  # noqa E402
from mantid.simpleapi import LoadNexusProcessed
from tempfile import NamedTemporaryFile


def _create_2d_data():
    """This creates the IQazimuthal data from ``Raw_Data_Anisotropic_v2_new`` the pages are

    * "Anisotropic Data - Qy vs Qx" contains the intensities with the Q values labeled in grey
    * "Anisotropic Data - EB Qy vs Qx" contains the uncertainties
    * "Anisotropic Data - Q" contains the Q-magnitude
    * "Anisotropic Data - Phi" contains the azimuthal angle
    """
    intensities = np.array(
        [
            [200, 190, 150, 190, 220, 210, 230, 250, 210, 190, 180],
            [180, 150, 220, 230, 230, 290, 190, 300, 180, 280, 310],
            [200, 200, 300, 210, 430, 380, 280, 290, 400, 380, 330],
            [400, 550, 420, 200, 220, 320, 320, 240, 700, 600, 350],
            [500, 600, 1100, 1500, 200, 180, 220, 1500, 1000, 700, 400],
            [600, 700, 1500, 2200, 3100, 0, 3000, 2300, 1300, 800, 500],
            [500, 600, 1200, 1500, 400, 250, 240, 1500, 1000, 700, 400],
            [400, 550, 230, 280, 380, 240, 200, 220, 700, 600, 350],
            [250, 320, 220, 220, 250, 340, 340, 290, 220, 220, 300],
            [200, 220, 180, 200, 300, 330, 230, 300, 180, 280, 320],
            [150, 190, 220, 180, 190, 280, 220, 290, 220, 220, 150],
        ]
    )
    errors = np.array(
        [
            [14.1, 13.8, 12.2, 13.8, 14.8, 14.5, 15.2, 15.8, 14.5, 13.8, 13.4],
            [13.4, 12.2, 14.8, 15.2, 15.2, 17.0, 13.8, 17.3, 13.4, 16.7, 17.6],
            [14.1, 14.1, 17.3, 14.5, 20.7, 19.5, 16.7, 17.0, 20.0, 19.5, 18.2],
            [20.0, 23.5, 20.5, 14.1, 14.8, 17.9, 17.9, 15.5, 26.5, 24.5, 18.7],
            [22.4, 24.5, 33.2, 38.7, 14.1, 13.4, 14.8, 38.7, 31.6, 26.5, 20.0],
            [24.5, 26.5, 38.7, 46.9, 55.7, 01.0, 54.8, 48.0, 36.1, 28.3, 22.4],
            [22.4, 24.5, 34.6, 38.7, 20.0, 15.8, 15.5, 38.7, 31.6, 26.5, 20.0],
            [20.0, 23.5, 15.2, 16.7, 19.5, 15.5, 14.1, 14.8, 26.5, 24.5, 18.7],
            [15.8, 17.9, 14.8, 14.8, 15.8, 18.4, 18.4, 17.0, 14.8, 14.8, 17.3],
            [14.1, 14.8, 13.4, 14.1, 17.3, 18.2, 15.2, 17.3, 13.4, 16.7, 17.9],
            [12.2, 13.8, 14.8, 13.4, 13.8, 16.7, 14.8, 17.0, 14.8, 14.8, 12.2],
        ]
    )

    # IQazimuthal's constructor is corrected in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/merge_requests/951
    # Thus the test input data shall be tranposed accordingly
    intensities = intensities.T
    errors = errors.T

    data2d = IQazimuthal(
        intensity=intensities,
        error=errors,
        qx=np.linspace(-5.0, 5.0, 11, dtype=float),
        qy=np.linspace(5.0, -5.0, 11, dtype=float),
    )
    assert data2d.intensity.shape == (11, 11)
    return data2d


def _create_2d_histogram_data():
    """This creates the parallel arrays of the binned data from ``Raw_Data_Anisotropic_v2_new`` the pages are

    * "Anisotropic Data - Q vs Phi" which contains the intensity of each Q column and the rows are the
      azimuthal angle values labeled with bin boundary. Both of these express bin centers.
    * "Anisotropic Data - EB Q vs Phi" which contains the uncertainties
    """
    # numbers taken from the spreadsheet
    intensity = np.array(
        [
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
            [3100, 2200, 1500, 700, 600, np.nan, np.nan],
        ],
        dtype=float,
    )

    error = np.array(
        [
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
            [55.7, 46.9, 38.7, 26.5, 24.5, np.nan, np.nan],
        ],
        dtype=float,
    )
    assert intensity.shape == error.shape

    # put together binning parameters to reproduce spreadsheet
    q_min = 0.5
    q_max = 7.5
    q_delta = 1.0
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)
    # verify the first and last values
    assert q_bins[0] == 0.5, "q_bins[0]"
    assert q_bins[-1] == 7.5, "q_bins[-1]"

    # create azimuthal angles as bin centers
    azimuthal_delta = 5.0
    azimuthal_max = 540.0
    azimuthal_bins = np.arange(
        start=0.0,
        stop=azimuthal_max + azimuthal_delta,
        step=azimuthal_delta,
        dtype=float,
    )
    # verify the first and last values
    assert azimuthal_bins[0] == 0.0, "azimuthal_bins[0]"
    assert azimuthal_bins[-1] == 540.0, "azimuthal_bins[-1]"
    azimuthal_bins += 2.5

    assert intensity.shape == (len(azimuthal_bins), len(q_bins) - 1)

    azimuthal_rings = []
    for i in range(intensity.shape[1]):
        azimuthal_rings.append(I1DAnnular(intensity=intensity.T[i], error=error.T[i], phi=azimuthal_bins))

    return q_bins, azimuthal_rings


def test_calc_qmod_and_azimuthal(cleanfile):
    """Test the conversion of data into 2d arrays of qmod and azimuthal angle. The results are
    checked against "Anisotropic Data -Q" and "Anisotropic Data - Phi"
    """
    data2d = _create_2d_data()

    # convert to q and azimuthal
    intensity, error, qmod, delta_qmod, azimuthal = _toQmodAndAzimuthal(data2d)
    assert qmod.shape == intensity.shape
    assert delta_qmod is None
    assert azimuthal.shape == intensity.shape

    from drtsans.plots.api import plot_IQazimuthal

    plot_IQazimuthal(data2d, "input_wedge.png", backend="mpl")
    cleanfile("input_wedge.png")

    # numbers taken from the spreadsheet
    q_exp = np.array(
        [
            [7.07, 6.40, 5.83, 5.39, 5.10, 5.00, 5.10, 5.39, 5.83, 6.40, 7.07],
            [6.40, 5.66, 5.00, 4.47, 4.12, 4.00, 4.12, 4.47, 5.00, 5.66, 6.40],
            [5.83, 5.00, 4.24, 3.61, 3.16, 3.00, 3.16, 3.61, 4.24, 5.00, 5.83],
            [5.39, 4.47, 3.61, 2.83, 2.24, 2.00, 2.24, 2.83, 3.61, 4.47, 5.39],
            [5.10, 4.12, 3.16, 2.24, 1.41, 1.00, 1.41, 2.24, 3.16, 4.12, 5.10],
            [5.00, 4.00, 3.00, 2.00, 1.00, 0.00, 1.00, 2.00, 3.00, 4.00, 5.00],
            [5.10, 4.12, 3.16, 2.24, 1.41, 1.00, 1.41, 2.24, 3.16, 4.12, 5.10],
            [5.39, 4.47, 3.61, 2.83, 2.24, 2.00, 2.24, 2.83, 3.61, 4.47, 5.39],
            [5.83, 5.00, 4.24, 3.61, 3.16, 3.00, 3.16, 3.61, 4.24, 5.00, 5.83],
            [6.40, 5.66, 5.00, 4.47, 4.12, 4.00, 4.12, 4.47, 5.00, 5.66, 6.40],
            [7.07, 6.40, 5.83, 5.39, 5.10, 5.00, 5.10, 5.39, 5.83, 6.40, 7.07],
        ],
        dtype=float,
    )

    azimuthal_exp = np.array(
        [
            [135, 129, 121, 112, 101, 90, 79, 68, 59, 51, 45],
            [141, 135, 127, 117, 104, 90, 76, 63, 53, 45, 39],
            [149, 143, 135, 124, 108, 90, 72, 56, 45, 37, 31],
            [158, 153, 146, 135, 117, 90, 63, 45, 34, 27, 22],
            [169, 166, 162, 153, 135, 90, 45, 27, 18, 14, 11],
            [180, 180, 180, 180, 180, 0, 0, 0, 0, 0, 0],
            [191, 194, 198, 207, 225, 270, 315, 333, 342, 346, 349],
            [202, 207, 214, 225, 243, 270, 297, 315, 326, 333, 338],
            [211, 217, 225, 236, 252, 270, 288, 304, 315, 323, 329],
            [219, 225, 233, 243, 256, 270, 284, 297, 307, 315, 321],
            [225, 231, 239, 248, 259, 270, 281, 292, 301, 309, 315],
        ],
        dtype=float,
    )

    np.testing.assert_allclose(qmod, q_exp.ravel(), atol=0.005)
    # IQazimuthal's constructor is corrected in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/merge_requests/951
    # Thus the expected test result shall be tranposed accordingly
    np.testing.assert_allclose(azimuthal, azimuthal_exp.T.ravel(), atol=0.5)


def test_bin_into_q_and_azimuthal():
    '''Test binning into Q and azimuthal matches the results from "Anisotropic Data - Q vs Phi"'''
    # get the test data
    data2d = _create_2d_data()
    q_exp, azimuthal_rings_exp = _create_2d_histogram_data()

    # parameters for azimuthal
    azimuthal_delta = 5.0

    # parameters for q
    q_min = 0.5
    q_max = 7.5
    q_delta = 1.0

    # get the histogrammed data
    q, azimuthal_rings = _binInQAndAzimuthal(
        data2d,
        q_min=q_min,
        q_max=q_max,
        q_delta=q_delta,
        azimuthal_delta=azimuthal_delta,
    )

    # verify the q-binning
    assert q.min() == q_min == q_exp.min()
    assert q.max() == q_max == q_exp.max()  # using bin boundaries

    for spectrum, spectrum_exp in zip(azimuthal_rings, azimuthal_rings_exp):
        assert spectrum.intensity.shape == spectrum_exp.intensity.shape
        np.testing.assert_allclose(spectrum.phi, spectrum_exp.phi, atol=0.05, equal_nan=True)

        np.testing.assert_allclose(spectrum.intensity, spectrum_exp.intensity, atol=0.05, equal_nan=True)
        np.testing.assert_allclose(spectrum.error, spectrum_exp.error, atol=0.05, equal_nan=True)


def test_fitting():
    """Test that the fitting generates reasonable results for fitting the peaks"""
    q, azimuthal_rings = _create_2d_histogram_data()
    # this calling forces there to be two found peaks
    selection = _fitQAndAzimuthal(
        azimuthal_rings,
        q,
        criteria=AcceptanceCriteria.from_settings(10.0, 0.2, 0.5, azimuthal_delta=5.0, maxchisq=1000.0),
        azimuthal_start=110.0,
        peak_search_window_size_factor=0.6,
    )
    assert selection.succeeded
    center_list, fwhm_list = selection.centers, selection.fwhms

    # Q-rings 2-5 are seed rings; ring 1 is rejected because its first peak has H/sigma_H of about 1
    assert center_list[0] == pytest.approx(180.0, abs=3.0)
    assert center_list[1] == pytest.approx(364.0, abs=4.0)
    assert fwhm_list[0] == pytest.approx(46.9, abs=2.0)
    assert fwhm_list[1] == pytest.approx(51.6, abs=2.0)
    # the two highest Q-rings have too few points to be fitted
    statuses = [ring.status for ring in selection.ring_fits]
    assert statuses == [RingStatus.REJECTED] + [RingStatus.SEED] * 4 + [RingStatus.FIT_FAILED] * 2
    assert all(ring.reasons for ring in selection.ring_fits if ring.status is not RingStatus.SEED)


def test_estimate_peak_parameters_empty_window():
    """An empty peak search window raises a RuntimeError with a clear message"""
    azimuthal = np.arange(0.0, 150.0, 5.0)  # too short to leave a window 90 deg inside both edges
    with pytest.raises(RuntimeError, match="no data in the peak search window"):
        _estimatePeakParameters(np.ones_like(azimuthal), azimuthal, azimuthal_start=110.0, window_half_width=108.0)


def _anisotropic_ring(contrast, background=1000.0, seed=42):
    """Azimuthal profile with two Gaussian peaks at 180 and 360 deg of height ``contrast`` x background"""
    phi = np.arange(2.5, 545.0, 5.0)
    peaks = sum(np.exp(-0.5 * (((phi - center + 180.0) % 360.0 - 180.0) / 15.0) ** 2) for center in (180.0, 360.0))
    intensity = background * (1.0 + contrast * peaks) + np.random.default_rng(seed).normal(0.0, 1.0, phi.size)
    return I1DAnnular(intensity=intensity, error=np.ones_like(phi), phi=phi)


def test_fit_low_contrast_ring():
    """A low-contrast ring is fitted, then rejected for its contrast only"""
    selection = _fitQAndAzimuthal(
        [_anisotropic_ring(contrast=0.1)],
        np.array([0.5, 1.5]),
        criteria=AcceptanceCriteria.from_settings(10.0, 0.2, 0.5, azimuthal_delta=5.0, maxchisq=1000.0),
        azimuthal_start=110.0,
        peak_search_window_size_factor=0.6,
    )
    (ring,) = selection.ring_fits
    assert ring.is_fitted and ring.status is RingStatus.REJECTED
    assert ring.reasons and all("H/A0" in reason for reason in ring.reasons if not reason.startswith("rescue"))
    assert selection.failure is SelectionFailure.NO_SEED_CANDIDATE


def test_get_wedge_selection_isotropic(tmp_path):
    """Nearly isotropic data raise AutoWedgeError, after the diagnostic files are written"""
    qx = qy = np.linspace(-5.0, 5.0, 41)
    qx_grid, qy_grid = np.meshgrid(qx, qy, indexing="ij")
    phi = np.arctan2(qy_grid, qx_grid)
    intensity = 1000.0 * (1.0 + 0.05 * np.cos(phi) ** 2) + np.random.default_rng(7).normal(0.0, 1.0, phi.shape)
    data2d = IQazimuthal(intensity=intensity, error=np.ones_like(intensity), qx=qx, qy=qy)

    with pytest.raises(AutoWedgeError, match="autoWedgeContrastMin") as error:
        getWedgeSelection(data2d, q_min=1.5, q_delta=1.0, q_max=4.5, azimuthal_delta=5.0, debug_dir=str(tmp_path))

    assert error.value.selection.failure is SelectionFailure.NO_SEED_CANDIDATE
    assert (tmp_path / "auto_wedge_rings.png").exists()
    assert (tmp_path / "auto_wedge_fit.h5").exists()


def test_fitting_symmetric():
    """Test that the fitting generates reasonable results for fitting the peaks"""
    q, azimuthal_rings = _create_2d_histogram_data()
    # this calling forces there to be two found peaks
    selection = _fitQAndAzimuthal(
        azimuthal_rings,
        q,
        criteria=AcceptanceCriteria.from_settings(10.0, 0.2, 0.5, azimuthal_delta=5.0, maxchisq=1000.0),
        azimuthal_start=110.0,
        peak_search_window_size_factor=0.6,
        auto_symmetric_wedges=True,
        auto_wedge_phi_min=2.5,
        auto_wedge_phi_max=217.5,
    )
    assert selection.succeeded
    center_list, fwhm_list = selection.centers, selection.fwhms

    assert center_list[0] == pytest.approx(180.0, abs=5.0)
    assert center_list[1] == pytest.approx(360.0, abs=5.0)
    assert fwhm_list[0] == pytest.approx(fwhm_list[1], abs=2.0)


def test_fitting_symmetric_negative_phi_min():
    """Test that a negative autoWedgePhiMin in symmetric mode gives the same peaks as the range shifted by 360 deg"""
    q, azimuthal_rings = _create_2d_histogram_data()
    criteria = AcceptanceCriteria.from_settings(10.0, 0.2, 0.5, azimuthal_delta=5.0, maxchisq=1000.0)
    selections = [
        _fitQAndAzimuthal(
            azimuthal_rings,
            q,
            criteria=criteria,
            azimuthal_start=110.0,
            peak_search_window_size_factor=0.6,
            auto_symmetric_wedges=True,
            auto_wedge_phi_min=phi_min,
            auto_wedge_phi_max=phi_min + 180.0,
        )  # fmt: skip
        for phi_min in (-90.0, 270.0)
    ]
    assert all(selection.succeeded for selection in selections)
    assert selections[0].centers == pytest.approx(selections[1].centers)
    assert selections[0].fwhms == pytest.approx(selections[1].fwhms)


def _wide_peaks_data() -> IQazimuthal:
    """Azimuthal profile 1 + 0.3 cos^2(phi) on a flat Q grid, with peaks of FWHM about 97 deg at 0 and 180 deg"""
    qx = qy = np.linspace(-5.0, 5.0, 41)
    qx_grid, qy_grid = np.meshgrid(qx, qy, indexing="ij")
    phi = np.arctan2(qy_grid, qx_grid)
    intensity = 1000.0 * (1.0 + 0.3 * np.cos(phi) ** 2) + np.random.default_rng(7).normal(0.0, 1.0, phi.shape)
    return IQazimuthal(intensity=intensity, error=np.ones_like(intensity), qx=qx, qy=qy)


def test_get_wedge_selection_no_room(tmp_path):
    """Test that peaks too wide for the background wedges at any background width raise an error"""
    # with autoWedgePeakWidth = 0.5, the reduced width S / (FWHM1 + FWHM2) - 0.5 is below the peak width
    with pytest.raises(AutoWedgeError, match="leave no room for background wedges") as error:
        getWedgeSelection(
            _wide_peaks_data(), q_min=1.5, q_delta=1.0, q_max=4.5, azimuthal_delta=5.0, peak_width=0.5,
            debug_dir=str(tmp_path),
        )  # fmt: skip
    assert error.value.selection.succeeded


def test_get_wedge_selection_wide_peaks(tmp_path, monkeypatch):
    """Peaks too wide for the requested background width give valid wedges and a warning"""
    data2d = _wide_peaks_data()
    logger = MagicMock()
    monkeypatch.setattr(drtsans.auto_wedge, "logger", logger)

    peak_wedges, background_wedges = getWedgeSelection(
        data2d, q_min=1.5, q_delta=1.0, q_max=4.5, azimuthal_delta=5.0, debug_dir=str(tmp_path)
    )

    (warning,), _ = logger.warning.call_args
    assert warning.startswith("autoWedgeBackgroundWidth reduced from 1.5 to")
    # the background wedge on the side of the smaller peak separation is as wide as the mean peak wedge,
    # the other one is wider by twice the deviation of the separation from 180 deg
    widths = [valid_wedge(*wedge) for wedge in peak_wedges + background_wedges]
    widths = [sum(high - low for low, high in pieces) for pieces in widths]
    assert min(widths[2:]) == pytest.approx(0.5 * (widths[0] + widths[1]), abs=0.5)
    centers = [0.5 * (low + high) for low, high in peak_wedges]
    separation = abs(centers[1] - centers[0]) % 360.0
    separation = min(separation, 360.0 - separation)
    assert max(widths[2:]) - min(widths[2:]) == pytest.approx(360.0 - 2.0 * separation, abs=0.5)


def test_get_wedge_selection_diagnostics(tmp_path):
    """The HDF5 file holds the status and reasons of every Q-ring, and the figure has one panel per ring"""
    getWedgeSelection(
        _create_2d_data(), q_min=0.5, q_max=7.5, q_delta=1.0, azimuthal_delta=5.0, debug_dir=str(tmp_path)
    )

    with h5py.File(tmp_path / "auto_wedge_fit.h5", "r") as h5:
        rings = {name: h5[group][name] for group in ("good fit", "no fit") for name in h5[group] if "ring" in name}
        statuses = [rings[f"ring {index}"]["status"][0].decode() for index in range(len(rings))]
        reasons = rings["ring 0"]["reasons"][0].decode().split("\n")
    assert statuses == ["rejected"] + ["seed"] * 4 + ["fit failed"] * 2
    assert reasons == ["peak 1: H/sigma_H = 1.02 < 10", "rescue: peak 1: H/sigma_H = 1.02 < 5"]

    # 7 Q-rings give 2 rows of panels
    width, height = RING_PLOT_PANEL_SIZE
    dpi = plt.rcParams["figure.dpi"]
    image = imread(tmp_path / "auto_wedge_rings.png")
    assert image.shape[:2] == (round(2 * height * dpi), round(RING_PLOT_COLUMNS * width * dpi))


def test_integration(tmp_path):
    """Test the full workflow of the algorithm"""
    data2d = _create_2d_data()
    # parameters for azimuthal
    azimuthal_delta = 5.0

    # parameters for q
    q_min = 0.5
    q_max = 7.5
    q_delta = 1.0

    # run the function
    wedges = getWedgeSelection(
        data2d,
        q_min=q_min,
        q_max=q_max,
        q_delta=q_delta,
        azimuthal_delta=azimuthal_delta,
        debug_dir=str(tmp_path),
    )

    # from drtsans.plots.api import plot_IQazimuthal
    # plot_IQazimuthal(data2d, 'wedge_int_test_debug.png', backend='mpl')
    # assert False, f'Wedges: {wedges}'

    # tests
    assert len(wedges) == 2
    for wedge in wedges:
        for min_val, max_val in wedge:
            print(min_val, "<", max_val)
            assert -90.0 < min_val < 270.0
            assert -90.0 < max_val < 270.0
    peak_wedge, back_wedge = wedges

    # Expected values follow from the Q-ring selection: rings 2-5 are seed rings, ring 1 is rejected because its
    # first peak has H/sigma_H of about 1. They differ from the spreadsheet, which used rings 1-4.
    # verify peaks
    # first peak
    assert 0.5 * (peak_wedge[0][0] + peak_wedge[0][1]) == pytest.approx(180.0, abs=3.0), "First peak center is at 180."
    assert peak_wedge[0][0] == pytest.approx(171.0, abs=0.5)
    assert peak_wedge[0][1] == pytest.approx(195.0, abs=0.5)
    # second peak
    assert 0.5 * (peak_wedge[1][0] + peak_wedge[1][1]) == pytest.approx(4.0, abs=1.0), "Second peak center is at 4."
    assert peak_wedge[1][0] == pytest.approx(-8.9, abs=0.5)
    assert peak_wedge[1][1] == pytest.approx(16.9, abs=0.5)

    # verify background
    # first background - the extra 360 is to get around the circle
    assert 0.5 * (back_wedge[0][0] + back_wedge[0][1] + 360) == pytest.approx(269.9, abs=1.2), (
        "First background center is at 270."
    )
    assert back_wedge[0][0] == pytest.approx(253.1, abs=0.5)
    assert back_wedge[0][1] == pytest.approx(-73.4, abs=0.5)
    # second background
    assert 0.5 * (back_wedge[1][0] + back_wedge[1][1]) == pytest.approx(97.0, abs=5.0), (
        "Second background center is at 97."
    )
    assert back_wedge[1][0] == pytest.approx(81.4, abs=0.5)
    assert back_wedge[1][1] == pytest.approx(112.6, abs=0.5)


def _assert_wedges(wedges, expected, tolerance=0.5):
    """Compare (peak, background) wedges with expected angle pairs, in degrees."""
    for found, wanted in zip(wedges, expected, strict=True):
        assert np.asarray(found) == pytest.approx(np.asarray(wanted), abs=tolerance)


@pytest.mark.datarepo
def test_real_data_eqsans(datarepo_dir, tmp_path):
    """Test the wedge selection on EQSANS data with a narrow streak on an incoherent background."""
    # the former mean/min ratio check rejected every Q-ring of these data
    data_file = os.path.join(datarepo_dir.eqsans, "auto_wedge", "EQSANS_auto_wedge_iqazimuthal.h5")
    iq2d = load_iq2d_from_h5(data_file)

    wedges = getWedgeSelection(
        iq2d, q_min=0.08, q_delta=0.05, q_max=0.8, azimuthal_delta=1.0, peak_width=0.25, debug_dir=str(tmp_path)
    )
    expected_peaks = [(176.6, 185.3), (-3.4, 5.9)]
    expected_backgrounds = [(207.1, -26.7), (29.2, 154.8)]
    _assert_wedges(wedges, (expected_peaks, expected_backgrounds))


@pytest.mark.datarepo
def test_real_data_eqsans_isotropic_q_range(datarepo_dir, tmp_path):
    """Test that the wedge selection fails with no seed ring on the isotropic Q range of the EQSANS data."""
    data_file = os.path.join(datarepo_dir.eqsans, "auto_wedge", "EQSANS_auto_wedge_iqazimuthal.h5")
    iq2d = load_iq2d_from_h5(data_file)

    with pytest.raises(AutoWedgeError, match="no Q-ring passes the seed criteria") as error:
        getWedgeSelection(
            iq2d, q_min=0.35, q_delta=0.05, q_max=0.8, azimuthal_delta=1.0, peak_width=0.25, debug_dir=str(tmp_path)
        )
    assert error.value.selection.failure is SelectionFailure.NO_SEED_CANDIDATE


@pytest.mark.datarepo
def test_real_data_eqsans_legacy_settings(datarepo_dir, tmp_path, monkeypatch):
    """Test that a background width equal to the peak width and a low significance give wedges and warnings."""
    data_file = os.path.join(datarepo_dir.eqsans, "auto_wedge", "EQSANS_auto_wedge_iqazimuthal.h5")
    iq2d = load_iq2d_from_h5(data_file)
    logger = MagicMock()
    monkeypatch.setattr(drtsans.auto_wedge, "logger", logger)

    peak_wedges, background_wedges = getWedgeSelection(
        iq2d, q_min=0.08, q_delta=0.05, q_max=0.8, azimuthal_delta=1.0, peak_width=0.25, background_width=0.25,
        signal_to_noise_min=2.0, debug_dir=str(tmp_path),
    )  # fmt: skip

    warnings = [args[0] for args, _ in logger.warning.call_args_list]
    assert any(warning.startswith("autoWedgeSignalToNoiseMin = 2 is the minimum significance") for warning in warnings)
    assert any("the background wedges overlap the peak wedges" in warning for warning in warnings)
    # adjacent wedges: each background wedge starts where a peak wedge ends
    assert background_wedges[0][0] == pytest.approx(peak_wedges[0][1])


@pytest.mark.datarepo
def test_real_data_biosans(datarepo_dir, tmp_path):
    """Test the wedge selection and the wedge binning on BIOSANS data."""
    MSamp_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5532_mBSub.h5")
    MBuff_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5562_mBSub.h5")

    ws_ms = LoadNexusProcessed(Filename=MSamp_fn, OutputWorkspace="sample", LoadHistory=False)
    ws_mb = LoadNexusProcessed(Filename=MBuff_fn, OutputWorkspace="Main_buffer", LoadHistory=False)
    ws_ms -= ws_mb  # subtract the buffer
    ws_mb.delete()

    # convert to I(qx,qy)
    q1d_data = biosans.convert_to_q(ws_ms, mode="scalar")
    q2d_data = biosans.convert_to_q(ws_ms, mode="azimuthal")
    ws_ms.delete()

    # calculate the wedge angles to use
    wedge_angles = getWedgeSelection(
        data2d=q2d_data,
        q_min=0.00,
        q_delta=0.001,
        q_max=0.02,
        azimuthal_delta=0.5,
        peak_width=0.25,
        background_width=1.5,
        debug_dir=str(tmp_path),
    )
    expected_peaks = [(157.3, 169.2), (-23.0, -10.7)]
    expected_backgrounds = [(198.9, -53.6), (19.9, 127.6)]
    _assert_wedges(wedge_angles, (expected_peaks, expected_backgrounds))

    # every wedge selects data, including the background wedges that wrap around 0 deg
    for wedge in wedge_angles[0] + wedge_angles[1]:
        pieces = [select_i_of_q_by_wedge(q2d_data, low, high) for low, high in valid_wedge(*wedge)]
        assert sum(piece.intensity.size for piece in pieces) > 0

    # test bin_all
    nbins = 100.0
    iq2d_rebinned, iq1d_rebinned = bin_all(
        q2d_data,
        q1d_data,
        nxbins=nbins,
        nybins=nbins,
        n1dbins=nbins,
        bin1d_type="wedge",
        wedges=wedge_angles,
        symmetric_wedges=False,
        error_weighted=False,
    )
    assert len(iq1d_rebinned) == 2, "Expect exactly 2 output 1d spectra"

    # rebin the data onto a regular grid for plotting
    linear_x_bins = determine_1d_linear_bins(q2d_data.qx.min(), q2d_data.qx.max(), nbins)
    linear_y_bins = determine_1d_linear_bins(q2d_data.qy.min(), q2d_data.qy.max(), nbins)
    q2d_data = bin_intensity_into_q2d(q2d_data, linear_x_bins, linear_y_bins, BinningMethod.NOWEIGHT)

    # save an image
    filename = NamedTemporaryFile(delete=False, prefix="CG3_127_5532_Iqxqy", suffix=".png").name
    plot_IQazimuthal(
        q2d_data,
        filename,
        backend="mpl",
        wedges=wedge_angles,
        symmetric_wedges=False,
        imshow_kwargs={"norm": LogNorm()},
    )
    print("saved image to", filename)
    plt.close()
    # verify the plot was created and remove the file
    assert os.path.exists(filename), '"{}" does not exist'.format(filename)
    os.remove(filename)


@pytest.mark.datarepo
def test_real_data_biosans_autosymmetric_wedges(datarepo_dir, tmp_path):
    """Test the symmetric wedge selection and the wedge binning on BIOSANS data."""
    MSamp_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5532_mBSub.h5")
    MBuff_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5562_mBSub.h5")

    ws_ms = LoadNexusProcessed(Filename=MSamp_fn, OutputWorkspace="sample", LoadHistory=False)
    ws_mb = LoadNexusProcessed(Filename=MBuff_fn, OutputWorkspace="Main_buffer", LoadHistory=False)
    ws_ms -= ws_mb  # subtract the buffer
    ws_mb.delete()

    # convert to I(qx,qy)
    q1d_data = biosans.convert_to_q(ws_ms, mode="scalar")
    q2d_data = biosans.convert_to_q(ws_ms, mode="azimuthal")
    ws_ms.delete()

    # calculate the wedge angles to use
    wedge_angles = getWedgeSelection(
        data2d=q2d_data,
        q_min=0.00,
        q_delta=0.001,
        q_max=0.02,
        azimuthal_delta=0.5,
        peak_width=0.25,
        background_width=1.5,
        auto_symmetric_wedges=True,
        debug_dir=str(tmp_path),
    )
    expected_peaks = [(-22.2, -11.5), (157.8, 168.5)]
    expected_backgrounds = [(15.4, 131.0), (195.4, -49.0)]
    _assert_wedges(wedge_angles, (expected_peaks, expected_backgrounds))

    # every wedge selects data, including the background wedges that wrap around 0 deg
    for wedge in wedge_angles[0] + wedge_angles[1]:
        pieces = [select_i_of_q_by_wedge(q2d_data, low, high) for low, high in valid_wedge(*wedge)]
        assert sum(piece.intensity.size for piece in pieces) > 0

    # test bin_all
    nbins = 100.0
    iq2d_rebinned, iq1d_rebinned = bin_all(
        q2d_data,
        q1d_data,
        nxbins=nbins,
        nybins=nbins,
        n1dbins=nbins,
        bin1d_type="wedge",
        wedges=wedge_angles,
        symmetric_wedges=False,
        error_weighted=False,
    )
    assert len(iq1d_rebinned) == 2, "Expect exactly 2 output 1d spectra"


@pytest.mark.datarepo
def test_real_data_biosans_manual(datarepo_dir):
    """Test asymmetric manual wedge binning on BIOSANS data"""
    MSamp_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5532_mBSub.h5")
    MBuff_fn = os.path.join(datarepo_dir.biosans, "CG3_127_5562_mBSub.h5")

    ws_ms = LoadNexusProcessed(Filename=MSamp_fn, OutputWorkspace="sample", LoadHistory=False)
    ws_mb = LoadNexusProcessed(Filename=MBuff_fn, OutputWorkspace="Main_buffer", LoadHistory=False)
    ws_ms -= ws_mb  # subtract the buffer
    ws_mb.delete()

    # convert to I(qx,qy)
    q1d_data = biosans.convert_to_q(ws_ms, mode="scalar")
    q2d_data = biosans.convert_to_q(ws_ms, mode="azimuthal")
    ws_ms.delete()

    # calculate the wedge angles to use
    wedge_angles = [(60.0, 120)]

    # test bin_all
    nbins = 100.0
    iq2d_rebinned, iq1d_rebinned = bin_all(
        q2d_data,
        q1d_data,
        nxbins=nbins,
        nybins=nbins,
        n1dbins=nbins,
        bin1d_type="wedge",
        wedges=wedge_angles,
        symmetric_wedges=False,
        error_weighted=False,
    )
    assert len(iq1d_rebinned) == 1, "Expect exactly 1 output 1d spectra for manual wedge"


if __name__ == "__main__":
    pytest.main([__file__])
