# Test drtsans.tof.eqsans.incoherence_correction_2d
import pytest
from drtsans.dataobjects import IQazimuthal
import numpy as np


def generate_test_data():
    # Generate test data
    # Once we have data anyway
    intensity_vec = np.array([np.nan])

    error_vec = np.sqrt(intensity_vec)

    vec_qx = np.array([np.nan])
    vec_qy = np.array([np.nan])

    wavelength_vec = np.array([np.nan])

    i_of_q = IQazimuthal(
        intensity=intensity_vec,
        error=error_vec,
        qx=vec_qx,
        qy=vec_qy,
        wavelength=wavelength_vec
    )

    return i_of_q
