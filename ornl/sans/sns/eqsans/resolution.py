from __future__ import (absolute_import, division, print_function)
import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.resolution import dq2_geometry, dq2_gravity
from ornl.sans.sns.eqsans import geometry as eqsans_geometry
from ornl.sans import geometry


def q_resolution_per_pixel(ws):
    r"""
    Compute q resolution for each pixel, in each wavelength bin.

    The resolution can be computed by giving a binned workspace to this function:

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
    L1 = geometry.sample_source_distance(ws, units='m')

    # TODO: the following would be ideal but the geometry implementation
    # needs logic to distinguish between time series and simple floats.
    # L2 = geometry.sample_detector_distance(ws,
    #                                        log_key='sample_detector_distance',
    #                                        units='m')
    sl = SampleLogs(ws)
    L2 = 1. / 1000. * sl.find_log_with_units('sample_detector_distance', 'mm')
    R1 = 1. / 2000. * eqsans_geometry.source_aperture(ws).diameter
    R2 = 1. / 2000. * eqsans_geometry.sample_aperture_diameter(ws)

    wl_bounds = ws.extractX()
    wl = (wl_bounds[:, 1:] + wl_bounds[:, :-1]) / 2.0
    dwl = wl_bounds[:, 1:] - wl_bounds[:, :-1]

    spec_info = ws.spectrumInfo()

    # Sanity check
    if not spec_info.size() == wl.shape[0]:
        raise RuntimeError("X size mismatch: %s %s" % \
                           (spec_info.size(), wl.shape[0]))

    theta = np.zeros_like(wl)
    qx = np.zeros_like(wl)
    qy = np.zeros_like(wl)
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            theta[i] = spec_info.twoTheta(i) * np.ones(wl.shape[1])
            _x, _y, _z = spec_info.position(i)
            phi = np.arctan2(_y, _x)
            _q = 4.0 * np.pi * np.sin(spec_info.twoTheta(i) / 2.0) / wl[i]
            qx[i] = np.cos(phi) * _q
            qy[i] = np.sin(phi) * _q

    dqx = np.sqrt(dqx2_eqsans(qx, L1, L2, R1, R2, wl, dwl, theta))
    dqy = np.sqrt(dqx2_eqsans(qy, L1, L2, R1, R2, wl, dwl, theta))
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
    if wl > 2.0:
        return 0.0148 * wl**3 - 0.5233 * wl**2 + 6.4797 * wl + 231.99
    else:
        return 392.31 * wl**6 - 3169.3 * wl**5 + 10445 * wl**4 \
            - 17872 * wl**3 + 16509 * wl**2 - 7448.4 * wl + 1280.5


def dqx2_eqsans(qx, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.011):
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
    return dq2_geo + np.fabs(qx) * (dwl / wl)**2 / 12.0


def dqy2_eqsans(qy, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.007):
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
    return dq2_geo + dq2_grav + np.fabs(qy) * (dwl / wl)**2 / 12.0
