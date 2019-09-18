import numpy as np

from drtsans.samplelogs import SampleLogs
from drtsans.momentum_transfer import dq2_geometry, dq2_gravity
from drtsans import geometry as sans_geometry

"""
HFIR Momentum transfer and resolution calculation
"""


def calculate_q_dq(ws, pixel_sizes=None):
    r"""
    Compute q resolution for each pixel, in each wavelength bin.

    The resolution can be computed by giving a binned
    workspace to this function:

    qx, qy, dqx, dqy = calculate_q_dq(ws_2d)

    The returned numpy arrays are of the same dimensions
    as the input array.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace
    pixel_sizes: None (use default) or dictionary of pixel size

    Returns
    ------
    2D arrays for Q, Qx, dQx, Qy, dQy, which are of the same dimension as the data
    """
    sl = SampleLogs(ws)
    L1 = sans_geometry.source_sample_distance(ws, unit='m',
                                              log_key='source-sample-distance')
    L2 = sans_geometry.sample_detector_distance(ws, unit='m')
    R1 = 0.5 * sl.find_log_with_units('source-aperture-diameter', 'mm') * 0.001
    R2 = 0.5 * sl.find_log_with_units('sample-aperture-diameter', 'mm') * 0.001
    wl = sl.find_log_with_units('wavelength')
    dwl = sl.find_log_with_units('wavelength-spread')

    # FIXME - Remove after testing
    print('Inputs: ')
    print("r1", R1)
    print("r2", R2)
    print("wl", wl)
    print("dwl", dwl)

    spec_info = ws.spectrumInfo()

    twotheta = np.zeros(spec_info.size())
    phi = np.zeros(spec_info.size())
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            twotheta[i] = spec_info.twoTheta(i)
            # assumes sample is at zero and orientation is McStas style
            _x, _y, _ = spec_info.position(i)
            phi[i] = np.arctan2(_y, _x)
        else:
            twotheta[i] = np.nan

    _q = 4.0 * np.pi * np.sin(0.5 * twotheta) / wl

    # convert things that are masked to zero
    _q[np.isnan(twotheta)] = 0.
    twotheta[np.isnan(twotheta)] = 0.  # do this one last

    qx = np.cos(phi) * _q
    qy = np.sin(phi) * _q
    del _q, phi

    print('Qx: {}'.format(qx))
    print('Qy: {}'.format(qy))
    print('L1: {}'.format(L1))
    print('L2: {}'.format(L2))
    print('R1: {}'.format(R1))
    print('R2: {}'.format(R2))
    print('WL: {}'.format(wl))
    print('dW: {}'.format(dwl))
    print('2T: {}'.format(twotheta))

    dqx = np.sqrt(_dqx2(qx, L1, L2, R1, R2, wl, dwl, twotheta))
    dqy = np.sqrt(_dqy2(qy, L1, L2, R1, R2, wl, dwl, twotheta))

    return qx, qy, dqx, dqy


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
        dq = np.sqrt(_dqx2(q_mid, L1, L2, R1, R2, wl, dwl))
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
            dqx[i] = np.sqrt(_dqx2(qx_mid, L1, L2, R1,
                                   R2, wl, dwl, theta))
            dqy[i] = np.sqrt(_dqy2(qy_mid[i], L1, L2,
                                   R1, R2, wl, dwl, theta))
        return dqx, dqy


def _dqx2(qx, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.0055):
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
    return dq2_geo + qx**2 * (dwl / wl)**2 / 6.0


def _dqy2(qy, L1, L2, R1, R2, wl, dwl, theta=None, pixel_size=0.0043):
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