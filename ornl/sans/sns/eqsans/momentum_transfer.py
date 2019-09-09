from __future__ import (absolute_import, division, print_function)
import numpy as np

from ornl.sans.momentum_transfer import dq2_geometry, dq2_gravity
from ornl.sans import geometry as sans_geometry
from ornl.sans.sns.eqsans import geometry as eqsans_geometry


def q_resolution_per_pixel(ws):
    """
    Compute q resolution for each pixel, in each wavelength bin.

    The resolution can be computed by giving a binned
    workspace to this function:

    qx, qy, dqx, dqy = q_resolution_per_pixel(ws_2d)

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
    L1 = sans_geometry.source_sample_distance(ws, unit='m',
                                              log_key='source-sample-distance')
    L2 = sans_geometry.sample_detector_distance(ws, unit='m')
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

    twotheta = np.zeros_like(wl)
    phi = np.zeros_like(wl)
    s2p = np.zeros_like(wl)
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            s2p[i] = spec_info.l2(i)
            twotheta[i] = spec_info.twoTheta(i)
            phi[i] = spec_info.azimuthal(i)
        else:
            twotheta[i] = np.nan

    mask = np.isnan(twotheta)
    _q = 4.0 * np.pi * np.sin(0.5 * twotheta) / wl
    _q[mask] = 0.
    twotheta[mask] = 0.  # do this one last
    del mask

    qx = np.cos(phi) * _q
    qy = np.sin(phi) * _q
    del _q, phi

    dtof = _moderator_time_error(wl)
    theta = 0.5 * twotheta
    dqx = np.sqrt(_dqx2(qx, L1, L2, R1, R2, wl, dwl, theta, s2p, dtof=dtof))
    dqy = np.sqrt(_dqy2(qy, L1, L2, R1, R2, wl, dwl, theta, s2p, dtof=dtof))
    return qx, qy, dqx, dqy


def _moderator_time_error(wl):
    """
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

    mask = wl > 2.0
    time_error[mask] = 0.0148 * wl[mask]**3 \
        - 0.5233 * wl[mask]**2 \
        + 6.4797 * wl[mask] + 231.99

    mask = wl <= 2.0
    time_error[mask] = 392.31 * wl[mask]**6 \
        - 3169.3 * wl[mask]**5 \
        + 10445 * wl[mask]**4 \
        - 17872 * wl[mask]**3 \
        + 16509 * wl[mask]**2 \
        - 7448.4 * wl[mask] + 1280.5

    del mask
    return time_error


def _dqx2(qx, L1, L2, R1, R2, wl, dwl, theta, s2p, pixel_size=0.0055, dtof=0.):
    """
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
    dq_tof_term = (3.9560 * dtof / 1000.0 / wl / (L1 + s2p))**2
    return dq2_geo + np.fabs(qx) * (dq_tof_term + (dwl / wl)**2) / 12.0


def _dqy2(qy, L1, L2, R1, R2, wl, dwl, theta, s2p, pixel_size=0.0043, dtof=0.):
    """
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
    dq_tof_term = (3.9560 * dtof / 1000.0 / wl / (L1 + s2p))**2
    dq2_grav = dq2_gravity(L1, L2, wl, dwl, theta)
    dq2 = dq2_geo + dq2_grav
    dq2 += np.fabs(qy) * (dq_tof_term + (dwl / wl)**2) / 12.0
    return dq2
