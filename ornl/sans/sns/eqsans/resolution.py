from __future__ import (absolute_import, division, print_function)
import numpy as np

from ornl.sans.resolution import dq2_geometry, dq2_gravity
from ornl.sans import geometry as sans_geometry
from ornl.sans.sns.eqsans import geometry as eqsans_geometry


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
    L1 = sans_geometry.source_sample_distance(ws, units='m',
                                              log_key='source-sample-distance')
    kwargs = dict(units='m', log_key='sample-detector-distance')
    L2 = sans_geometry.sample_detector_distance(ws, **kwargs)
    R1 = 0.5 * eqsans_geometry.source_aperture_diameter(ws, unit='m')
    R2 = 0.5 * eqsans_geometry.sample_aperture_diameter(ws, unit='m')

    wl_bounds = ws.extractX()
    wl = (wl_bounds[:, 1:] + wl_bounds[:, :-1]) / 2.0
    dwl = wl_bounds[:, 1:] - wl_bounds[:, :-1]

    spec_info = ws.spectrumInfo()

    # Sanity check
    if not spec_info.size() == wl.shape[0]:
        raise RuntimeError("X size mismatch: %s %s" %
                           (spec_info.size(), wl.shape[0]))

    theta = np.zeros_like(wl)
    qx = np.zeros_like(wl)
    qy = np.zeros_like(wl)
    s2p = np.zeros_like(wl)
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            theta[i] = spec_info.twoTheta(i) * np.ones(wl.shape[1])
            s2p[i] = spec_info.l2(i) * np.ones(wl.shape[1])
            _x, _y, _z = spec_info.position(i)
            phi = np.arctan2(_y, _x)
            _q = 4.0 * np.pi * np.sin(spec_info.twoTheta(i) / 2.0) / wl[i]
            qx[i] = np.cos(phi) * _q
            qy[i] = np.sin(phi) * _q

    dqx = np.sqrt(dqx2_eqsans(qx, L1, L2, R1, R2, wl, dwl, theta, s2p))
    dqy = np.sqrt(dqx2_eqsans(qy, L1, L2, R1, R2, wl, dwl, theta, s2p))
    return dqx, dqy


def moderator_time_error(wl):
    r"""
    Relative Q uncertainty due to emission time jitter
    in the neutron moderator.

    Parameters
    ----------
    wl: float
        wavelength [Angstrom]

    Returns
    float
    """
    time_error = np.zeros_like(wl)
    time_error[wl > 2.0] = 0.0148 * wl[wl > 2.0]**3 \
        - 0.5233 * wl[wl > 2.0]**2 \
        + 6.4797 * wl[wl > 2.0] + 231.99
    time_error[wl <= 2.0] = 392.31 * wl[wl <= 2.0]**6 \
        - 3169.3 * wl[wl <= 2.0]**5 \
        + 10445 * wl[wl <= 2.0]**4 \
        - 17872 * wl[wl <= 2.0]**3 \
        + 16509 * wl[wl <= 2.0]**2 \
        - 7448.4 * wl[wl <= 2.0] + 1280.5
    return time_error


def dqx2_eqsans(qx, L1, L2, R1, R2, wl, dwl, theta, s2p, pixel_size=0.011):
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
    s2p: float
        sample-to-pixel (m)
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
    dtof = moderator_time_error(wl)
    dq_tof_term = (3.9560 * dtof / 1000.0 / wl / (L1 + s2p))**2
    return dq2_geo + np.fabs(qx) * (dq_tof_term + (dwl / wl)**2) / 12.0


def dqy2_eqsans(qy, L1, L2, R1, R2, wl, dwl, theta, s2p, pixel_size=0.007):
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
    s2p: float
        sample-to-pixel (m)
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
    dtof = moderator_time_error(wl)
    dq_tof_term = (3.9560 * dtof / 1000.0 / wl / (L1 + s2p))**2
    dq2_grav = dq2_gravity(L1, L2, wl, dwl, theta)
    dq2 = dq2_geo + dq2_grav
    dq2 += np.fabs(qy) * (dq_tof_term + (dwl / wl)**2) / 12.0
    return dq2
