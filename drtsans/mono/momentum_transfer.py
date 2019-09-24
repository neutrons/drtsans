import numpy as np

from drtsans.samplelogs import SampleLogs
from drtsans.momentum_transfer import dq2_geometry, dq2_gravity, InstrumentSetupParameters, calculate_momentum_transfer
from drtsans import geometry as sans_geometry

"""
HFIR Momentum transfer and resolution calculation
"""


def calculate_q_resolution(qx, qy,  wave_length, delta_wave_length, theta, instrument_setup_params):
    """ Calculate Q resolution
    Parameters
    ----------
    qx : float or float array
        Qx or float array
    qy : float
        Qy
    wave_length : float or float array
        wave length in A
    delta_wave_length : float or float array
        wave length in A
    theta : float or float array
        theta value
    instrument_setup_params : InstrumentSetupParameters
        collection of instrument setup related parameters including     l1, l2, r1, r2, pixel_size_x, pixel_size_y

    Returns
    -------
    (Float, Float) or (numpy.array, numpy.array)
        dQx, dQy
    """
    # # If theta is not supplied, compute it from qx
    # # This simplifies the calculation for I(Q) in 1D.
    # if theta is None:
    #     # FIXME - one of it must be wrong!
    #     thetax = 2.0 * np.arcsin(wave_length * np.fabs(qx) / 4.0 / np.pi)
    #     thetay = 2.0 * np.arcsin(wave_length * np.fabs(qy) / 4.0 / np.pi)
    #     assert thetax == thetay

    print('[DEBUG] WL = {} +/- {}\nTheta = {}\n{}'.format(wave_length, delta_wave_length,
                                                          theta, instrument_setup_params))

    # Calculate dQx
    dqx = np.sqrt(_dqx2(qx,
                        instrument_setup_params.l1,
                        instrument_setup_params.sample_det_center_distance,
                        instrument_setup_params.source_aperture_radius,
                        instrument_setup_params.sample_aperture_radius,
                        wave_length, delta_wave_length, theta,
                        instrument_setup_params.pixel_size_x))

    dqy = np.sqrt(_dqy2(qy, instrument_setup_params.l1,
                        instrument_setup_params.sample_det_center_distance,
                        instrument_setup_params.source_aperture_radius,
                        instrument_setup_params.sample_aperture_radius,
                        wave_length, delta_wave_length, theta,
                        instrument_setup_params.pixel_size_y))

    print('----> dQx = {}, dQy = {}'.format(dqx, dqy))

    return dqx, dqy


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
    ndarray, ndarray, ndarray, ndarry
        Qx, Qy, dQx, dQy
        2D arrays for Qx, dQx, Qy, dQy, which are of the same dimension as the data
    """
    # Check inputs and unit
    if ws is None:
        raise RuntimeError('Workspace cannot be None')
    elif ws.getAxis(0).getUnit().unitID() != 'Wavelength':
        raise RuntimeError('Input workspace {} for calculate Q resolution must be in unit Wavelength but not {}'
                           ''.format(ws, ws.getAxis(0).getUni().unitID()))
    elif not ws.isHistogramData():
        raise RuntimeError('Input workspace {} must be Point Data but not histogram data.'.format(ws))

    # Get instrument setup parameters to calculate Q resolution
    wl, dwl, setup_params = retrieve_instrument_setup(ws, pixel_sizes)

    # Calculate momentum_transfer Q
    q, qx, qy, two_theta, s2p = calculate_momentum_transfer(ws)

    # Calculate Q resolution
    dqx, dqy = calculate_q_resolution(qx, qy, wl, dwl, two_theta*0.5, setup_params)

    return qx, qy, dqx, dqy


def retrieve_instrument_setup(ws, pixel_sizes):
    """

    Parameters
    ----------
    ws
    pixel_sizes

    Returns
    -------
    Float, Float, InstrumentSetupParameters
        wave length, delta wave length, InstrumentSetupParameters instance (L1, L2, ...)

    """
    # Retrieve sample logs
    sl = SampleLogs(ws)

    # L1
    l1 = sans_geometry.source_sample_distance(ws, unit='m',
                                              log_key='source-sample-distance')
    # L2
    l2 = sans_geometry.sample_detector_distance(ws, unit='m')
    # R1: source aperture in meter
    r1 = 0.5 * sl.find_log_with_units('source-aperture-diameter', 'mm') * 0.001
    # R2: sample aperture in meter
    r2 = 0.5 * sl.find_log_with_units('sample-aperture-diameter', 'mm') * 0.001
    # Constant wave length and uncertainties are recorded as sample logs in HFIR SANS
    wl = sl.find_log_with_units('wavelength')
    dwl = sl.find_log_with_units('wavelength-spread')

    # Pixel sizes
    if pixel_sizes is None:
        # Retrieve from workspace but not easy
        det_shape = ws.getDetector(0).shape().getBoundingBox().width()  # 3 values
        size_x = det_shape[0]
        size_y = det_shape[1]
    else:
        # User specified, overriding values from intrument directly
        size_x = pixel_sizes['x']
        size_y = pixel_sizes['y']

    # Set up the parameter class
    setup_params = InstrumentSetupParameters(l1=l1,
                                             sample_det_center_dist=l2,
                                             source_aperture_radius=r1,
                                             sample_aperture_radius=r2,
                                             pixel_size_x=size_x,
                                             pixel_size_y=size_y)

    print('[DEBUG INFO] Instrument parameters:\nlambda = {} +/- {}\n{}'
          ''.format(wl, dwl, setup_params))

    return wl, dwl, setup_params


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


def _dqx2(qx, L1, L2, R1, R2, wl, dwl, theta, pixel_size=0.0055):
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
    dq2_geo = dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size)
    dq2_wl = qx**2 * (dwl / wl)**2 / 6.0

    if isinstance(qx, np.ndarray):
        print('[DEBUG....SPECIAL] Qx   : shape = {} value = \n{}'
              ''.format(qx.shape, qx))
        print('[DEBUG....SPECIAL] Dq2_geo: shape = {} value = \n{}'
              ''.format(dq2_geo.shape, dq2_geo))
        print('[DEBUG....SPECIAL] Dq2_wl:  shape = {} value = \n{}'
              ''.format(dq2_wl.shape, dq2_wl))

    return dq2_geo + dq2_wl


def _dqy2(qy, L1, L2, R1, R2, wl, dwl, theta, pixel_size=0.0043):
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
    dq2_geo = dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size)
    dq2_grav = dq2_gravity(L1, L2, wl, dwl, theta)

    return dq2_geo + dq2_grav + qy**2 * (dwl / wl)**2 / 6.0
