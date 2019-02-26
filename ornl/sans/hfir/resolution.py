from __future__ import (absolute_import, division, print_function)
import numpy as np

from ornl.sans.samplelogs import SampleLogs


def q_resolution(ws):
    r"""
    Compute q resolution for the given reduced workspace.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace

    Returns
    ------
    numpy array of the same dimension as the data
    """
    sl = SampleLogs(ws)

    # TODO: the following two lines don't work for HFIR
    # source_sample = geometry.source_sample_distance(ws, units='m')
    # sample_detector = geometry.sample_detector_distance(ws, units='m')

    L1 = 1. / 1000. * sl.find_log_with_units('source-sample-distance', 'mm')
    L2 = 1. / 1000. * sl.find_log_with_units('sample-detector-distance', 'mm')
    R1 = 1. / 2000. * sl.find_log_with_units('source-aperture-diameter', 'mm')
    R2 = 1. / 2000. * sl.find_log_with_units('sample-aperture-diameter', 'mm')
    wl = sl.find_log_with_units('wavelength', 'Angstrom')
    dwl = sl.find_log_with_units('wavelength-spread', 'Angstrom')

    q = ws.extractX()
    r = ws.extractY()
    if len(q) == 1:
        # We have a 1D I(q)
        q_mid = (q[0][1:] + q[0][:-1]) / 2.0
        dq = np.sqrt(dqy2_hfir(q_mid, L1, L2, R1, R2, wl, dwl))
        return dq
    else:
        # We have a 2D I(qx, qy)
        # TODO: make sure that the correct axis is Y
        dqx = np.zeros_like(r)
        dqy = np.zeros_like(r)
        _qx_values = np.asarray(ws.getAxis(0).extractValues())
        _qy_values = np.asarray(ws.getAxis(1).extractValues())
        qx_mid = (_qx_values[1:] + _qx_values[:-1]) / 2.0
        qy_mid = (_qy_values[1:] + _qy_values[:-1]) / 2.0
        for i in range(len(q)):
            q_length = np.sqrt(qy_mid[i]**2 + qx_mid**2)
            theta = 2.0 * np.arcsin(wl * q_length / 4.0 / np.pi)
            dqx[i] = np.sqrt(dqx2_hfir(qx_mid, L1, L2, R1,
                                       R2, wl, dwl, theta))
            dqy[i] = np.sqrt(dqy2_hfir(qy_mid[i], L1, L2,
                                       R1, R2, wl, dwl, theta))
        return dqx, dqy


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
    dq2 = 0.25 * (L2 / L1 * R1)**2 + 0.25 * ((L1 + L2) / L1 * R2)**2 \
        + pixel_size**2 / 12.0
    return dq2 * (2.0 * np.pi * np.cos(theta) \
                  * np.cos(2.0 * theta)**2 / wl / L2)**2


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
    h = 6.626e-34    # m^2 kg s^-1
    m_n = 1.675e-27  # kg
    g = 9.8          # m s^-2
    B = g * m_n**2 * L2 * (L1 + L2) / 2 / h**2
    dq2 = 2 / 3 * B**2 * wl**2 * dwl**2
    dq2 *= (2.0 * np.pi * np.cos(theta) * np.cos(2.0 * theta)**2 / wl / L2)**2
    # Converting from A^2 / m^4 to 1 / A^2
    return dq2 * 1.0e-40


def dqx2_hfir(qx, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.011):
    r"""
    Q resolution in the horizontal direction.

    Parameters
    ----------
    qx: float
        value of q_x (1/Angstrom)
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
    # If theta is not supplied, compute it from qx
    # This simplifies the calculation for I(Q) in 1D.
    if theta is None:
        theta = 2.0 * np.arcsin(wl * np.fabs(qx) / 4.0 / np.pi)
    dq2_geo = dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size)
    return dq2_geo + np.fabs(qx) * (dwl / wl)**2 / 6.0


def dqy2_hfir(qy, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.007):
    r"""
    Q resolution in vertical direction.

    Parameters
    ----------
    qy: float
        value of q_y (1/Angstrom)
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
    # If theta is not supplied, compute it from qx
    # This simplifies the calculation for I(Q) in 1D.
    if theta is None:
        theta = 2.0 * np.arcsin(wl * np.fabs(qy) / 4.0 / np.pi)
    dq2_geo = dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size)
    dq2_grav = dq2_gravity(L1, L2, wl, dwl, theta)
    return dq2_geo + dq2_grav + np.fabs(qy) * (dwl / wl)**2 / 6.0
