import numpy as np

from drtsans.momentum_transfer import dq2_geometry, dq2_gravity
from drtsans import geometry as sans_geometry
from drtsans.tof.eqsans.geometry import source_aperture_diameter, sample_aperture_diameter, source_sample_distance
from drtsans.momentum_transfer import calculate_momentum_transfer, MomentumTransferResolutionParameters
from drtsans.momentum_transfer import G_MN2_OVER_H2


def calculate_q_dq(ws, pixel_sizes=None):
    """
    Calculate momentum transfer and momentum transfer resolution for each pixel in each wavelength bin

    The resolution can be computed by giving a binned
    workspace to this function:

    qx, qy, dqx, dqy = calculate_q_dq(ws_2d)

    The returned numpy arrays are of the same dimensions
    as the input array.

    Note:
    1. pixel sizes: the auto-retrieved value of detector pixel size comes from det.boundBox(), which may be slightly
                    larger than the real pixel size.  For example: 0.0042972564697265625 vs. 0.004296875 for height.
                    User can specify the pixel sizes to override.


    :param ws: MatrixWorkspace in unit of wave length.
               It must be histogram data for bin step (point data does not have this information)
    :param pixel_sizes: 2-tuple for size x and size y of each pixel or None (default)
    :return: numpy array (shape = (n, m) n = ws.getNumberHistograms(), m = number of wave length bins) for Qx, Qy,
             dQx, dQy
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
    setup_params = retrieve_instrument_setup(ws, pixel_sizes)

    # From histogram to get wave length at center of each bin. The output is a 2D array for pixel number and
    # wave length bin
    # FIXME - This is only good for constant binning
    wavelength_bin_boundary_matrix = ws.extractX()
    wavelength_bin_center_matrix = 0.5 * (wavelength_bin_boundary_matrix[:, 1:] +
                                          wavelength_bin_boundary_matrix[:, :-1])
    wavelength_bin_step_matrix = wavelength_bin_boundary_matrix[:, 1:] - wavelength_bin_boundary_matrix[:, :-1]

    print('[DEBUG.........] wave length matrix: {}, {}, {}'.format(wavelength_bin_boundary_matrix.shape,
                                                                   wavelength_bin_center_matrix.shape,
                                                                   wavelength_bin_step_matrix.shape))

    # Calculate detector pixel information
    # TODO - late in the issue: move this part to calculate_momentum_transfer() to reduce loop
    vectors = calculate_pixel_positions(ws)
    pixel_2theta_vec, pixel_sample_distance_vec = vectors
    # pixel_phi_vec
    print('[DEBUG.........] s2p, 2theta: {}, {}'.format(pixel_sample_distance_vec .shape, pixel_2theta_vec.shape))

    # Calculate momentum transfer Q
    returns = calculate_momentum_transfer(ws)
    qx_matrix = returns[1]
    qy_matrix = returns[2]

    # Calculate neutron emission time according to wave length value
    tof_error_matrix = moderator_time_uncertainty(wavelength_bin_center_matrix)

    # Calculate dQx and dQy, both as 2D arrays
    dqx_matrix, dqy_matrix = calculate_q_resolution(qx_matrix, qy_matrix,
                                                    wavelength_bin_center_matrix, wavelength_bin_step_matrix,
                                                    0.5 * pixel_2theta_vec, pixel_2theta_vec,
                                                    pixel_sample_distance_vec,
                                                    tof_error_matrix, setup_params)

    return qx_matrix, qy_matrix, dqx_matrix, dqy_matrix


def q_resolution_per_pixel_to_mod(ws):
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
    L1 = source_sample_distance(ws, unit='m', log_key='source-sample-distance')
    L2 = sans_geometry.sample_detector_distance(ws, unit='m')
    R1 = 0.5 * source_aperture_diameter(ws, unit='m')
    R2 = 0.5 * sample_aperture_diameter(ws, unit='m')

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

    dtof = moderator_time_uncertainty(wl)
    theta = 0.5 * twotheta
    dqx = np.sqrt(_dqx2(qx, L1, L2, R1, R2, wl, dwl, theta, s2p, dtof=dtof))
    dqy = np.sqrt(_dqy2(qy, L1, L2, R1, R2, wl, dwl, theta, s2p, dtof=dtof))
    return qx, qy, dqx, dqy


def moderator_time_uncertainty(wl):
    """ Relative Q uncertainty due to emission time jitter
    in the neutron moderator.
    :param wl: float (or ndarray) wavelength [Angstrom]
    :return: float or ndarray (same shape to wave_length_array) of emission error time
    """
    # init output array to zeros
    time_error = np.zeros_like(wl)

    # formula for lambda>2
    mask = wl > 2.0
    time_error[mask] = 0.0148 * wl[mask]**3 \
        - 0.5233 * wl[mask]**2 \
        + 6.4797 * wl[mask] + 231.99

    # formula for lambda<=2
    mask = wl <= 2.0
    time_error[mask] = 392.31 * wl[mask]**6 \
        - 3169.3 * wl[mask]**5 \
        + 10445 * wl[mask]**4 \
        - 17872 * wl[mask]**3 \
        + 16509 * wl[mask]**2 \
        - 7448.4 * wl[mask] + 1280.5

    # clean up memory
    del mask

    # convert from zero-dimensional nparray to float
    if isinstance(wl, float):
        time_error = float(time_error)

    return time_error


def retrieve_instrument_setup(ws, pixel_sizes=None):
    """ Get instrument parameter including L1, L2, source aperture diameter and sample aperture radius
    :param ws:
    :param pixel_sizes: dictionary for pixel sizes
    :return: MomentumTransferResolutionParameters instance
    """
    # Retrieve L1 and L2 from instrument geometry
    l1 = source_sample_distance(ws, unit='m', search_logs=False)
    l2 = sans_geometry.sample_detector_distance(ws, unit='m',
                                                search_logs=False)
    r1 = 0.5 * source_aperture_diameter(ws, unit='m')
    r2 = 0.5 * sample_aperture_diameter(ws, unit='m')

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
    setup_params = MomentumTransferResolutionParameters(l1=l1,
                                                        sample_det_center_dist=l2,
                                                        source_aperture_radius=r1,
                                                        sample_aperture_radius=r2,
                                                        pixel_size_x=size_x,
                                                        pixel_size_y=size_y)

    return setup_params


def calculate_pixel_positions(ws, num_spec=None):
    """
    Get pixel 2theta, theta, pixel-sample-distance
    :param ws:
    :param num_spec: number of spectra
    :return: 2 vectors of pixels' (2theta, distance to sample)
    """
    # Get each spectrum
    spec_info = ws.spectrumInfo()

    # Sanity check
    if num_spec is None:
        num_spec = spec_info.size()
    elif num_spec != spec_info.size():
        raise RuntimeError('Number of spectra do not match!')

    # Initialize arrays
    pixel_2theta_vec = np.zeros(shape=(num_spec, 1), dtype='float')
    pixel_sample_distance_vec = np.zeros_like(pixel_2theta_vec)

    # Calculate theta/2theta and will be very slow!
    for i in range(num_spec):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            # sample pixel distance
            pixel_sample_distance_vec[i, 0] = spec_info.l2(i)
            # 2theta
            pixel_2theta_vec[i, 0] = spec_info.twoTheta(i)  # unit = radians
        else:
            # no detector or monitor (not likely)
            pixel_2theta_vec[i, 0] = np.nan
            pixel_sample_distance_vec[i, 0] = np.nan
        # END-IF-ELSE
    # END-FOR

    return pixel_2theta_vec, pixel_sample_distance_vec


def calculate_q_resolution(qx, qy, wave_length, delta_wave_length, theta, two_theta, sample_pixel_distance,
                           tof_error, q_resolution_params):
    """ Atomic function to calculate Q resolution for EQ-SANS
    :param qx: Qx
    :param qy: Qy
    :param wave_length: neutron wave length (bin center) in Angstrom
    :param delta_wave_length: neutron wave length (bin size) in Angstrom
    :param theta: half neutron diffraction angle (half of 2theta) (unit: rad)
    :param two_theta: neutron diffraction angle (unit: rad)
    :param sample_pixel_distance: distance from sample to pixel center (meter)
    :param tof_error: neutron emission uncertainty
    :param q_resolution_params: MomentumTransferResolutionParameters parameters
    :return:
    """
    assert isinstance(q_resolution_params, MomentumTransferResolutionParameters)
    print('Input Q resolution parameter:\n{}'.format(q_resolution_params))
    print('Q = {}, {}'.format(qx, qy))
    print('Wavelength = {}, Delta Wavelength = {}'.format(wave_length, delta_wave_length))
    print('Theta = {}, 2Theta = {}'.format(theta, two_theta))
    print('TOF uncertainty = {}'.format(tof_error))
    # Get setup value (for better)

    l1 = q_resolution_params.l1
    l2 = q_resolution_params.sample_det_center_distance
    r1 = q_resolution_params.source_aperture_radius
    r2 = q_resolution_params.sample_aperture_radius
    dx = q_resolution_params.pixel_size_x
    dy = q_resolution_params.pixel_size_y

    # Geometry resolution/uncertainty:
    apertures_const = (l2*r1*0.5/l1)**2 + (0.5 * r2 * (l1 + l2) / l1)**2
    const_x = dx**2 / 12.
    const_y = dy**2 / 12.

    # factor 1 with ....
    factor1 = (2.*np.pi*np.cos(theta)*(np.cos(two_theta)**2)/(wave_length*l2))**2
    print('[EQ Res] Factor1 = {}'.format(factor1))
    print('[EQ Res] Resolution X = {}'.format(apertures_const + const_x))
    # FIXME - Good so far

    # gravity part of Y ....
    b_factor = 0.5 * G_MN2_OVER_H2 * (l1 + l2) * l2
    gravity_y = 2./3. * (b_factor * wave_length * delta_wave_length * 1.E-20)**2  # convert to meter

    print('[EQ Res] B = {}'.format(b_factor))
    print('[EQ Res] Resolution Y (gravity) = {}'.format(gravity_y))

    # FIXME - Good so far

    # wave length resolution/uncertainties
    wl_part = (delta_wave_length/wave_length)**2
    emission_part = (3.9650*0.001*tof_error / (wave_length*(l1+sample_pixel_distance)))**2
    print('[DEBUG INFO] wave length uncertainty square = {}; emission time uncertainty square = {}'
          ''.format(wl_part, emission_part))

    print('[DEBUG INFO] factor 1 shape = {}'.format(factor1.shape))

    # Q(x) resolution
    qx_geom_resolution = factor1 * (apertures_const + const_x)
    qx_wave_resolution = qx**2 * (wl_part + emission_part) / 12.
    dqx = qx_wave_resolution + qx_wave_resolution

    print('Resolution X = {}'.format(apertures_const + const_x))

    # Q(y) resolution
    qy_geom_resolution = factor1 * (apertures_const + const_y + gravity_y)
    qy_wave_resolution = qy**2 * (wl_part + emission_part) / 12.
    dqy = qy_geom_resolution + qy_wave_resolution

    print('[DEBUG INFO] Geom dQx = {}, Wave dQx = {}'.format(qx_geom_resolution, qx_wave_resolution))
    print('[DEBUG INFO] Geom dQy = {}, Wave dQy = {}'.format(qy_geom_resolution, qy_wave_resolution))

    return dqx, dqy


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
    return dq2_geo + qx**2 * (dq_tof_term + (dwl / wl)**2) / 12.0


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
    dq2 += qy**2 * (dq_tof_term + (dwl / wl)**2) / 12.0
    return dq2
