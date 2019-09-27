import numpy as np
from scipy import constants
from drtsans.momentum_transfer import dq2_geometry, dq2_gravity
from drtsans import geometry as sans_geometry
from drtsans.tof.eqsans.geometry import source_aperture_diameter, sample_aperture_diameter, source_sample_distance
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
from mantid.simpleapi import LoadEmptyInstrument, AddTimeSeriesLog, Rebin, ConvertUnits
from drtsans.tof.eqsans.momentum_transfer import calculate_q_dq, calculate_pixel_positions,\
    retrieve_instrument_setup
from drtsans.tof.eqsans.momentum_transfer import calculate_q_resolution, InstrumentSetupParameters,\
    moderator_time_uncertainty
from drtsans.tof.eqsans import load_events


# This test implements issue #187 to verify methods working for q_resolution_per_pixel
# q_resolution_per_pixel() is tested in test_resolution
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Wei-Ren Chen


def assert_delta(exp_value, gold_value, delta, value_name):
    """
    Assert float number equal with delta
    :param exp_value:
    :param gold_value:
    :param delta:
    :param value_name:
    :return:
    """
    if abs(exp_value - gold_value) > delta:
        raise AssertionError('{}: Expected = {}, Testing = {} Different over {}'
                             ''.format(value_name, gold_value, exp_value, delta))

    return


def assert_equal(exp_value, gold_value, value_name):
    """
    Assert float number equal with delta
    :param exp_value:
    :param gold_value:
    :param value_name:
    :return:
    """
    if exp_value != gold_value:
        raise AssertionError('{}: Expected = {}, Testing = {}'
                             ''.format(value_name, gold_value, exp_value))

    return


def _set_generic_workspace(ws):
    # TODO : Get a list of wave length: valid range: [ 2.46179456  6.71611136] A

    # Range of TOF
    wave_length_range = np.array([2.5, 4.5])  # A
    intensity = np.array([20.]*16) + np.random.randn(16,)
    init_delta_intensity = np.random.randn(intensity.shape[0], )

    # assume that the TOF is already frame corrected
    for i in range(16):
        ws.dataX(i)[:] = wave_length_range  # microseconds
        ws.dataY(i)[0] = intensity[i]
        ws.dataE(i)[0] = init_delta_intensity[i]
    # #### ABOVE THIS POINT WILL BE A TEST FIXTURE

    # Set apertures' logs
    ws = _set_apertures_log(ws)

    return ws


def _set_apertures_log(input_ws):
    """ Add the sample logs to
    :param input_ws:
    :return:
    """
    # Add logs for source aperture (3) and Sample aperture (1)
    log_names = ['vBeamSlit', 'vBeamSlit2', 'vBeamSlit3', 'beamslit4']
    log_values = [5., 5., 5., 10.]  # all in unit mm
    input_ws_name = input_ws.name()
    for ilog in range(len(log_names)):
        log_name = log_names[ilog]
        log_value = log_values[ilog]
        AddTimeSeriesLog(input_ws, Name=log_name, Time="2010-01-01T00:00:00", Value=log_value)

    if input_ws is None:
        raise RuntimeError('Input {}/{} is not None; but returned is'.format(input_ws, input_ws_name))

    return input_ws


def check_pixels_position(ws):
    """ check pixels' position on the generic workspace
    :param ws:
    :return:
    """
    if ws.getNumberHistograms() != 16:
        raise AssertionError('Workspace {} shall have 16 spectra but not {}'
                             ''.format(ws, ws.getNumberHistograms()))

    # Golden position value
    gold_pos_dict = {0: np.array([0.009, -0.006, 1.25]),  # lower_left_pos
                     6: np.array([0.003,  0.002, 1.24]),  # middle (left and upper)
                     12: np.array([-0.009, -0.006, 1.25]),  # lower_right_pos
                     15:  np.array([-0.009, 0.006, 1.25])   # upper_right_pos
                     }

    # Check 3 corners
    match = True
    for iws in [0,  12, 15]:
        if np.sqrt(np.sum((ws.getDetector(iws).getPos() - gold_pos_dict[iws])**2)) > 1.E-10:
            match = False
            break

    # Assertion
    if not match:
        pixel_pos_str = ''
        for iws in range(16):
            pixel_pos_str += '{}: {}\n'.format(iws, ws.getDetector(iws).getPos())
        raise AssertionError(pixel_pos_str)

    # Get pixel positions
    outputs = calculate_pixel_positions(ws, ws.getNumberHistograms())
    pixel_2theta_vec, pixel_sample_distance_vec = outputs

    # Verify 2theta (diffraction angle) in radians: 2 corners with 2theta and pixel sample distance
    # 2theta = np.atan(np.sqrt(pos_0**2 + pos_1**2) / pos_2**2))  pos_2 == L2 = 1.25
    assert_delta(pixel_2theta_vec[0],  0.008653107083873311, 1.E-10, '2theta @ lower left corner (0)')
    assert_delta(pixel_2theta_vec[6],  0.002884433020894695, 1.E-10, '2theta @ middle lower upper (6)')
    assert_delta(pixel_2theta_vec[15], 0.008653107083873311, 1.E-10, '2theta @ upper right corner (15)')

    # Verify pixel-L2: 2 corners with 2theta and pixel sample distance
    # np.sqrt(np.sum(det_pos - [0])**2)
    assert_delta(pixel_sample_distance_vec[0], 1.2500467991239368, 1.E-10, 'Sample-to-lower left corner (0)')
    assert_delta(pixel_sample_distance_vec[6], 1.2500051999891841, 1.E-10, 'Sample-to-left-upper middle (6)')
    assert_delta(pixel_sample_distance_vec[15], 1.2500467991239368, 1.E-10, 'Sample-to-upper right corner (15)')

    return


# TODO FIXME - Disabled as EQSANS_92160 cannot be found!
def not_test_info_retrieve_real_nexus():
    """
    Test instrument parameters and meta data retrieval from real EQSANS data
    by recognizing the special SANS IDF
    :return:
    """
    eqsans_ws = load_events('EQSANS_92160')  # MetaDataOnly=True

    # Convert unit and rebin
    eqsans_ws = ConvertUnits(InputWorkspace=eqsans_ws, Target='Wavelength')
    eqsans_ws = Rebin(eqsans_ws, Params='-0.001')

    # Get L1, L2, R1 and R2
    setup_dict = retrieve_instrument_setup(eqsans_ws, pixel_sizes=None)

    # {'R1': 0.01, 'R2': 0.005, 'L2': 1.3, 'L1': 14.122}
    assert_delta(setup_dict.l1, 14.122, 1.E-7, 'L1')
    assert_delta(setup_dict.sample_det_center_distance, 1.30, 1.E-7, 'L2')
    assert_delta(setup_dict.source_aperture_radius, 0.0075, 1.E-7, 'R1 (source aperture)')
    # 0.5 * (0.005 + 0.005 + 0.005)
    assert_delta(setup_dict.sample_aperture_radius, 0.005, 1.E-7, 'R2 (sample aperture')
    # 0.5 * 10 mm
    # '<ns1:radius val="0.0055"/>\n      <ns1:height val="0.004296875"/>\n
    #    </ns1:cylinder>\n
    assert_delta(setup_dict.pixel_size_x, 0.011, 1.E-7, 'Pixel X size')  # 0.006 meter
    assert_delta(setup_dict.pixel_size_y, 0.004296875, 1.E-6, 'Pixel Y size')
    # bound box on height (0.0042972564697265625 is slightly larger)

    return


def test_single_value_resolution():
    """ Test method calculate Q resolution mostly from Generic SANS
    :return:
    """
    # l1 = 14.122
    # l2 = 1.30
    # source_aperture = 0.0075  # source aperture
    # sample_aperture = 0.0050  # sample aperture
    # qx = 0.01553395  # lambda 3.5 corner pixel  (0, 1)
    # qy = 2.48543251e-05  # lambda 3.5  corner pixel (1, 1)
    # wave_length = 3.5
    # wl_resolution = 0.5
    # two_theta = 0.008653107083873311  # radian (corner pixel)
    # sample_pixel_distance = 1.2500467991239368  # radian (corner pixel)
    # emission_error = 248.893075  # wave length = 3.5 A

    """
    Q = -0.000593411755, -0.000767944624
    Wavelength = 6.0, Delta Wavelength = 0.15
    Theta = 0.00046338, 2Theta = 0.00092676
    Pixel size = 0.0055, 0.0043
    R1 = 0.02, R2 = 0.007
    L1 = 15, L2 = 15.5
    B = 0.00014811607134644243
    """

    l1 = 15.
    l2 = 15.5
    source_aperture = 0.02  # source aperture
    sample_aperture = 0.007  # sample aperture
    qx = -0.000593411755
    qy = -0.000767944624
    wave_length = 6.0
    wl_resolution = 0.15
    two_theta = 0.00092676  # radian (corner pixel)
    sample_pixel_distance = l2 + 0.1  # radian (corner pixel)
    emission_error = 0.  # wave length = 3.5 A

    params = InstrumentSetupParameters(l1=l1,
                                       sample_det_center_dist=l2,
                                       source_aperture_radius=source_aperture,
                                       sample_aperture_radius=sample_aperture,
                                       pixel_size_x=0.0055,
                                       pixel_size_y=0.0043)

    q_x_res, q_y_res = calculate_q_resolution(qx=qx, qy=qy, wave_length=wave_length, delta_wave_length=wl_resolution,
                                              two_theta=two_theta,
                                              sample_pixel_distance=sample_pixel_distance,
                                              tof_error=emission_error,
                                              instrument_setup_params=params)
    print('Unit Test: Backend dQx = {}'.format(q_x_res))

    # Calculte Q resolution by Weiren's algorithm
    golden_dqx, golden_dqy = sigma_neutron(wave_length, wl_resolution, qx, qy, 0.5*two_theta,
                                           l1, l2, source_aperture, sample_aperture, 0.0055, 0.0043,
                                           sample_pixel_distance, l1, emission_error)

    assert_delta(q_x_res, golden_dqx, 1E-12, 'Q_x resolution')
    assert_delta(q_y_res, golden_dqy, 1E-12, 'Q_y resolution')

    return


@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4,
                           'dx': 0.006, 'dy': 0.004, 'zc': 1.25,  # TODO - it is best to use close to real L1 and L2
                           'l1': -5.}],
                         indirect=True)
def skip_test_q_resolution_per_pixel(generic_IDF):
    """
    Create a generic (SANS) instrument and test for main()
    :param generic_IDF:
    :return:
    """
    # Generate a generic SANS instrument with a pixel of
    # the size and position specified in
    # sans-backend/documents/Master_document_022219.pdf
    with open(r'/tmp/GenericSANS_Definition.xml', 'w') as tmp:
        tmp.write(generic_IDF)
        tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS',
                             OutputWorkspace='test_uncertainty')
    ws.getAxis(0).setUnit('Wavelength')

    # Set sample logs and other testing information
    ws = _set_generic_workspace(ws)

    # Check pixel's positions
    check_pixels_position(ws)

    # Set uncertainties
    # qx_matrix, qy_matrix, dqx_matrix, dqy_matrix = calculate_q_dq(ws)
    calculate_q_dq(ws)

    # Qx and Qy shall be already tested in sans.momentum_transfer

    return


def sigma_neutron(wave_length, delta_wave_length, Qx, Qy, theta, L1, L2, R1, R2, x3, y3, s2p, m2s, sig_emission):
    """
    Q resolution calculation from Wei-ren and modified by wzz
    Parameters
    ----------
    wave_length: float
        neutron wavelength (A)
    delta_wave_length: float
        wavelength width (A)
    Qx: float
        momentum transfer in x direction
    Qy: float
         momentum transfer in y direction
    theta: float
        scattering angle (half of 2theta of the pixel)
    L1: float
        source to sample
    L2: float
        sample to center of detector
    R1: float
        source aperture (m)
    R2: float
        sample aperture (m)
    x3: float
        detector pixel dimensions along x-axis
    y3: float
        detector pixel dimensions along y-axis
    s2p: float
        sample to pixel distance (m), which is slightly different from L2
    m2s: float
        moderator to sample distance
    sig_emission: float
        neutron emission time (second?)

    Returns
    -------
    float, float
        dQx, dQy

    """
    # For EQ-SANS

    # Define constants
    h = 6.62607004e-34  # scipy.const.h
    mn = 1.674929e-27  # scipy.const.neutron_mass
    g = constants.g  # 9.81 n/s^2

    # Calculate B
    B = 0.5 * g * mn**2 * L2 * (L1+L2) / h**2
    B /= 10**20  # add a factor of 10^-2 to cancel out the A in the term using B
    # (dWL/WL)**2
    r = (delta_wave_length/wave_length)**2

    print('[UNIT TEST WeiRen] B = {}, r = {}'.format(B, r))

    # dQx
    sigma_x = (2. * np.pi * np.cos(theta) * np.cos(2.*theta)**2 / wave_length / L2)**2
    sigma_x = sigma_x * ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2*R2**2/4 + x3**2/12)  # geometry
    print('[UNIT TEST WeiRen] Geometry dQx = {}'.format(sigma_x))
    sigma_x = sigma_x + Qx**2 / 12 * (r + (3.9560*sig_emission)**2/(1000*wave_length*(s2p+m2s))**2)
    print('[UNIT TEST WeiRen] Wavelength dQx = {}'
          ''.format(Qx**2 / 12 * (r + (3.9560*sig_emission)**2/(1000*wave_length*(s2p+m2s))**2)))

    # dQy
    sigma_y = (2. * np.pi * np.cos(theta) * np.cos(2*theta)**2 / wave_length / L2)**2
    print('[UNIT TEST WeiRen] Qy/wl factor = {}'.format(sigma_y))
    print('[UNIT TEST WeiRen] Geometry factor = {}, {}, {} = {}'
          ''.format((L2/L1)**2*R1**2/4 + (1+L2/L1)**2*R2**2/4, y3 ** 2 / 12,
                    B ** 2 * wave_length ** 4 * 2 / 3 * r,
                    ((L2 / L1) ** 2 * R1 ** 2 / 4 + (1 + L2 / L1) ** 2 * R2 ** 2 / 4 +
                     y3 ** 2 / 12 + B ** 2 * wave_length ** 4 * 2 / 3 * r)))
    sigma_y = sigma_y * ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2*R2**2/4 + y3**2/12 + B**2*wave_length**4*2/3*r)  # geometry
    print('[UNIT TEST WeiRen] Geometry dQy = {}'.format(sigma_y))
    sigma_y = sigma_y + Qy**2 / 12 * (r + (3.9560*sig_emission)**2/(1000*wave_length*(s2p+m2s))**2)
    print('[UNIT TEST WeiRen] Wavelength dQy = {}'
          ''.format(Qy**2 / 12 * (r + (3.9560*sig_emission)**2/(1000*wave_length*(s2p+m2s))**2)))

    print('sigma_x = {:.9E} sigma_y = {:.9E} sigma = {:.9f}\n'
          ''.format(sigma_x, sigma_y, np.sqrt(sigma_x**2+sigma_y**2)))

    return sigma_x, sigma_y


# FIXME TODO - put this to 2 parts...
#              1. Qx and Qy calculation
#              2. Q resolution calculation
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
