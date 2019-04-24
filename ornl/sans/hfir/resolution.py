from __future__ import (absolute_import, division, print_function)
import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.resolution import dq2_geometry, dq2_gravity
from ornl.sans import geometry as sans_geometry

def q_resolution_per_pixel(ws):
    r"""
    Compute q resolution for each pixel, in each wavelength bin.

    The resolution can be computed by giving a binned
    workspace to this function:

    dqx, dqy = q_resolution_per_pixel(ws_2d)

    The returned numpy arrays are of the same dimensions
    as the input array.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace

    Returns
    ------
    numpy array of the same dimension as the data
    """
    sl = SampleLogs(ws)
    L1 = sans_geometry.source_sample_distance(ws, 'm')
    L2 = sans_geometry.sample_detector_distance(ws, 'm')
    R1 = 1. / 2000. * sl.find_log_with_units('source-aperture-diameter', 'mm')
    R2 = 1. / 2000. * sl.find_log_with_units('sample-aperture-diameter', 'mm')
    wl = sl.find_log_with_units('wavelength', 'Angstrom')
    dwl = sl.find_log_with_units('wavelength-spread', 'Angstrom')

    spec_info = ws.spectrumInfo()

    theta = np.zeros(spec_info.size())
    qx = np.zeros(spec_info.size())
    qy = np.zeros(spec_info.size())
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            theta[i] = spec_info.twoTheta(i)
            _x, _y, _z = spec_info.position(i)
            phi = np.arctan2(_y, _x)
            _q = 4.0 * np.pi * np.sin(spec_info.twoTheta(i) / 2.0) / wl
            qx[i] = np.cos(phi) * _q
            qy[i] = np.sin(phi) * _q

    dqx = np.sqrt(dqx2_hfir(qx, L1, L2, R1, R2, wl, dwl, theta))
    dqy = np.sqrt(dqx2_hfir(qy, L1, L2, R1, R2, wl, dwl, theta))
    return dqx, dqy


def q_resolution(ws):
    r"""
    Compute q resolution for the given reduced workspace.

    The resolution can be computed by giving either a 1D I(q) or
    2D I(qx, qy) to this function:

    dq = q_resolution(ws_1d)
    dqx, dqy = q_resolution(ws_2d)

    In both cases, the returned numpy arrays are of the same dimensions
    as the input array.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace

    Returns
    ------
    numpy array of the same dimension as the data
    """
    sl = SampleLogs(ws)

    L1 = sans_geometry.source_sample_distance(ws, 'm')
    L2 = sans_geometry.sample_detector_distance(ws, 'm')
    R1 = 1. / 2000. * sl.find_log_with_units('source-aperture-diameter', 'mm')
    R2 = 1. / 2000. * sl.find_log_with_units('sample-aperture-diameter', 'mm')
    wl = sl.find_log_with_units('wavelength', 'Angstrom')
    dwl = sl.find_log_with_units('wavelength-spread', 'Angstrom')

    q = ws.extractX()
    r = ws.extractY()
    if len(q) == 1:
        # We have a 1D I(q)
        q_mid = (q[0][1:] + q[0][:-1]) / 2.0
        dq = np.sqrt(dqx2_hfir(q_mid, L1, L2, R1, R2, wl, dwl))
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
