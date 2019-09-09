"""
    Resolution calculations common to all SANS
"""
from __future__ import absolute_import, division, print_function

import numpy as np
from scipy import constants

# derived constant where:
# h = 6.626e-34    # m^2 kg s^-1
# m_n = 1.675e-27  # kg
# g = 9.8          # m s^-2
_G_MN2_OVER_H2 = constants.G * np.square(constants.neutron_mass / constants.h)


def dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size=0.007):
    r"""
    Geometry term for both dqx and dqy.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace
    L1: float
        source-to-sample distance (m)
    L2: float
        sample-to-detector distance (m)
    R1: float
        source aperture radius (m)
    R2: float
        sample aperture radius (m)
    wl: float
        wavelength mid-point (Angstrom)
    theta: float
        scattering angle (rad)
    pixel_size: float
        dimension of the pixel (m)

    Returns
    ------
    float
    """
    dq2 = 0.25 * np.square(L2 / L1 * R1)  \
        + 0.25 * np.square((L1 + L2) / L1 * R2) + np.square(pixel_size) / 12.0
    return dq2 * np.square(2.0 * np.pi * np.cos(theta)
                           * np.square(np.cos(2.0 * theta)) / wl / L2)


def dq2_gravity(L1, L2, wl, dwl, theta):
    r"""
    Gravity term for dqy.

    Parameters
    ----------
    L1: float
        source-to-sample distance (m)
    L2: float
        sample-to-detector distance (m)
    wl: float
        wavelength mid-point (Angstrom)
    dwl: float
        wavelength-spread (Angstrom)
    theta: float
        scattering angle (rad)
    pixel_size: float
        dimension of the pixel (m)

    Returns
    ------
    float
    """
    B = _G_MN2_OVER_H2 * 0.5 * L2 * (L1 + L2)
    dq2 = 2. * np.square(B * wl * dwl) / 3.
    dq2 *= np.square(2.0 * np.pi * np.cos(theta)
                     * np.square(np.cos(2.0 * theta)) / wl / L2)
    # Converting from A^2 / m^4 to 1 / A^2
    return dq2 * 1.0e-40
