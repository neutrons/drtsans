import numpy as np
import pytest

from mantid.simpleapi import LoadHFIRSANS, AddSampleLog
from drtsans.momentum_transfer import InstrumentSetupParameters
from drtsans.mono.momentum_transfer import calculate_q_dq, calculate_q_resolution
from drtsans.samplelogs import SampleLogs
import scipy
import scipy.constants


# This implements Issue #168: calculate dQx and dQy
# dev - Wenduo Zhou <wzz@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Wei-Ren Chen


def sigma_neutron_weiren(wavelength, delta_lambda, Qx, Qy, theta, L1, L2, R1, R2, x3, y3):
    """
    Function given by Wei-Ren

    input: sigma_neutron_weiren(7, 0.7, 0.03, 0.03, 0.1, 10, 20, 5, 10, 0.03, 0.03, 20, 50)


    For GP-SANS and Bio-SANS
    lambda: netron wavelength in angstrom
    delta_lambda: wavelength width in angstrom
    Qx,Qy: momentum transfer in x and y direction, respectively, with unit (1/angstrom)
    theta: scattering angle in radian
    L1,L2: the flight path lengths whose units are meter.
    R1, R2: sample and source apertures, respectively, in meter.
    x3,y3: detector pixel dimensions in meter

    """
    h = 6.62607004e-34
    # h = scipy.constants.h
    mn = 1.674929e-27
    # mn = scipy.constants.neutron_mass
    g = scipy.constants.g  # 6.67408e-11
    B = 0.5*g*mn**2*L2*(L1+L2)/h**2
    B = B / 10**20
    r = (delta_lambda/wavelength)**2
    sigma_x = (2*np.pi*np.cos(theta)*np.cos(2*theta)**2 / wavelength/L2)**2
    sigma_x = sigma_x * ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2 * R2**2/4 + x3**2 / 12)
    sigma_x = sigma_x + Qx**2/6*r
    sigma_y = (2*np.pi*np.cos(theta)*np.cos(2*theta)**2/wavelength/L2)**2
    sigma_y = sigma_y * ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2 *
                         R2**2/4 + y3**2/12 + 2*B**2*wavelength**4
                         * r/3) + Qy**2/6*r
    sigma_x = np.sqrt(sigma_x)
    sigma_y = np.sqrt(sigma_y)

    return sigma_x, sigma_y


def sigma_neutron(wavelength, delta_lambda, Qx, Qy, theta, L1, L2, R1, R2, x3, y3):
    """
    Function given by Wei-Ren and rewritten by WZZ for detailed comparison

    input: sigma_neutron_weiren(7, 0.7, 0.03, 0.03, 0.1, 10, 20, 5, 10, 0.03, 0.03, 20, 50)


    For GP-SANS and Bio-SANS
    lambda: netron wavelength in angstrom
    delta_lambda: wavelength width in angstrom
    Qx,Qy: momentum transfer in x and y direction, respectively, with unit (1/angstrom)
    theta: scattering angle in radian
    L1,L2: the flight path lengths whose units are meter.
    R1, R2: sample and source apertures, respectively, in meter.
    x3,y3: detector pixel dimensions in meter

    """
    h = scipy.constants.h
    mn = scipy.constants.neutron_mass
    g = scipy.constants.g  # 6.67408e-11
    B = 0.5*g*mn**2*L2*(L1+L2)/h**2
    B = B / 10**20
    r = (delta_lambda/wavelength)**2

    print('Q = {}, {}'.format(Qx, Qy))
    # print('Wavelength = {}, Delta Wavelength = {}'.format(wavelength, delta_lambda))
    # print('Theta = {}, 2Theta = {}'.format(theta, 2*theta))
    # print('Pixel size = {}, {}'.format(x3, y3))
    # print('R1 = {}, R2 = {}'.format(R1, R2))
    # print('L1 = {}, L2 = {}'.format(L1, L2))
    # print('B = {}'.format(B))

    print('[UNIT TEST]WL = {} +/- {}, Theta = {}, L1 = {}, L2 = {}, R1 = {}, R2 = {}'
          ', Pixel = {}, {}'.format(wavelength, delta_lambda, theta, L1, L2,
                                    R1, R2, x3, y3))

    # Qy resolution
    sigma_x = (2*np.pi*np.cos(theta)*np.cos(2*theta)**2 / wavelength/L2)**2
    resolution_x = ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2 * R2**2/4 + x3**2 / 12)

    # print('X: factor1 = {}, resolution = {}'.format(sigma_x, resolution_x))
    sigma_x = sigma_x * resolution_x
    print('[UNIT TEST] dQ2_geo = {}, dQ2_Qx = {}'.format(sigma_x, Qx**2 * r / 6.))
    sigma_x = sigma_x + Qx**2 * r / 6.
    sigma_x = np.sqrt(sigma_x)

    # Qy resolution
    sigma_y = (2*np.pi*np.cos(theta)*np.cos(2*theta)**2/wavelength/L2)**2
    resolution_y = ((L2/L1)**2*R1**2/4 + (1+L2/L1)**2 * R2**2/4 + y3**2/12 + 2*B**2*wavelength**4 * r/3)
    # print('Y gravity part = {}'.format(2*B**2*wavelength**4 * r/3))
    # print('Y parr1 = {}'.format(sigma_y*resolution_y))
    sigma_y = sigma_y * resolution_y + Qy**2/6*r
    sigma_y = np.sqrt(sigma_y)
    # print('sigma_x = {:.2e}; sigma_y = {:.2e}; sigma = {:.2e}'.format(
    #     sigma_x, sigma_y, sigma_x+sigma_y))

    return sigma_x, sigma_y


def test_sigma_neutron():
    """ Test whether the rewritten sigma_neutron is same as Wei-ren's
    """
    R1 = 0.02
    R2 = 0.007

    x3 = 0.0055
    y3 = 0.0043
    L1 = 15
    L2 = 15.5

    wavelength = 6.0
    delta_lambda = 0.15

    Qx = -5.93411755e-04
    Qy = -7.67944624e-04

    # two_theta = 0.00092676    # radian (corner pixel)
    theta = 0.00092676 * 0.5  # radian (corner pixel)

    dqx, dqy = sigma_neutron_weiren(wavelength, delta_lambda, Qx, Qy, theta, L1, L2, R1, R2, x3, y3)

    dqx2, dqy2 = sigma_neutron(wavelength, delta_lambda, Qx, Qy, theta, L1, L2, R1, R2, x3, y3)

    assert dqx == pytest.approx(dqx2, 1.E-12)
    assert dqy == pytest.approx(dqy2, 1.E-12)

    # Calculate by drtsans.mono method
    setup = InstrumentSetupParameters(L1, L2, R1, R2, x3, y3)
    dqx3, dqy3 = calculate_q_resolution(Qx, Qy, wavelength, delta_lambda, theta, setup)

    # Check
    assert dqx == pytest.approx(dqx3, 1E-12)
    assert dqy == pytest.approx(dqy3, 1E-12)

    # assert 'Helping EQSANS Test' == 'Pass'

    return


@pytest.mark.parametrize('generic_workspace', [{
    'name': 'GPSANS',
    'Nx': 5, 'Ny': 5,
    'dx': 0.00425, 'dy': 0.0055,
    'xc': 0.0, 'yc': 0.0, 'zc': 15.5, 'l1': 15,
    'axis_values': [5.925, 6.075],
}], indirect=True)
def test_q_resolution_generic(generic_workspace):
    """
    Test Q resolution method against Wei-ren and Ricardo's early implementation
    Parameters
    ----------
    generic_workspace : Workspace instance
        A generic workspace with 5 x 5 instrument

    Returns
    -------

    """
    # Define constants
    wavelength = 6
    delta_lambda = 0.15
    R1 = 0.02  # source aperture radius
    R2 = 0.007  # sample aperture radius
    x3 = 0.00425  # pixel X size (meter)
    y3 = 0.00550  # pixel Y size (meter)
    L1 = 15  # meter
    L2 = 15.5  # meter (sample to detetor center distance)

    # Get workspace and add sample logs as wave length, wave length spread,
    ws = generic_workspace
    AddSampleLog(Workspace=ws, LogName='wavelength', LogText='{}'.format(wavelength), LogType='Number', LogUnit='A')
    AddSampleLog(Workspace=ws, LogName='wavelength-spread', LogText='{}'.format(delta_lambda),
                 LogType='Number', LogUnit='A')
    AddSampleLog(Workspace=ws, LogName='source-aperture-diameter', LogText='{}'.format(R1*2.*1000),
                 LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=ws, LogName='sample-aperture-diameter', LogText='{}'.format(R2*2.*1000),
                 LogType='Number', LogUnit='mm')

    # Calculate Q and dQ
    qx_arr, qy_arr, dqx_arr, dqy_arr = calculate_q_dq(ws)
    # Check dimension
    assert qx_arr.shape == (25, 1)
    assert dqy_arr.shape == (25, 1)

    # Calculate Q and 2theta (arrays)
    ver_qx_array, ver_qy_array, two_theta_array = calculate_momentum_transfer(ws, wavelength)

    # Check by another approach!
    spec_info = ws.spectrumInfo()
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            # compare Qx and Qy
            qx = qx_arr[i][0]
            qy = qy_arr[i][0]

            qx_check = ver_qx_array[i]
            qy_check = ver_qy_array[i]

            print(qx, qx_check, qx-qx_check)
            print(qy, qy_check, qy-qy_check)

            assert abs(qx - qx_check) < 1E-8, 'Spectrum {} Qx error: {} vs {} (test version)' \
                                              ''.format(i, qx, qx_check)
            assert abs(qy - qy_check) < 1E-8, 'Spectrum {} Qy error'.format(i)

            assert qx == pytest.approx(qx_check, 1.E-6)
            assert qy == pytest.approx(qy_check, 1.E-6)

            # calculate resolution
            # print('[UNIT TEST] Spectrum {}: WL = {} +/- {}, Theta = {}, L1 = {}, L2 = {}, R1 = {}, R2 = {}'
            #       ', Pixel = {}, {}'.format(i, wavelength, delta_lambda, 0.5*two_theta_array, L1, L2,
            #                                 R1, R2, x3, y3))
            sigma_x, sigma_y = sigma_neutron(wavelength, delta_lambda, qx, qy, 0.5*two_theta_array[i],
                                             L1, L2, R1, R2, x3, y3)
            assert abs(sigma_x - dqx_arr[i][0]) < 1E-08, 'Spectrum {} dQx error'.format(i)

            # assert sigma_x == pytest.approx(dqx_arr[i][0], 1.E-6)
            # assert sigma_y == pytest.approx(dqy_arr[i][0], 1.E-6)
        # END-IF
    # END-FOR

    return


def calculate_momentum_transfer(ws, wl):
    """Previous implementation for calculating Q by using Phi and Theta
    Parameters
    ----------
    ws : workspace instance
        Workspace where Q is calculated from
    wl:  float
        (constant) wave length

    Returns
    -------
    ndarray, ndarray, ndarray
        Qx, Qy, 2theta
    """
    spec_info = ws.spectrumInfo()
    twotheta = np.zeros(spec_info.size())
    phi = np.zeros(spec_info.size())
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            twotheta[i] = spec_info.twoTheta(i)
            phi_i = spec_info.azimuthal(i)
            phi[i] = phi_i
        else:
            twotheta[i] = np.nan
    # END-FOR

    print('Phi: shape = {}\n{}'.format(phi.shape, phi.reshape((25, 1))))

    _q = 4.0 * np.pi * np.sin(0.5 * twotheta) / wl

    # convert things that are masked to zero
    _q[np.isnan(twotheta)] = 0.
    twotheta[np.isnan(twotheta)] = 0.  # do this one last

    print('Q:  shape = {}\n{}'.format(_q.shape, _q.reshape((25, 1))))

    qx = np.cos(phi) * _q
    qy = np.sin(phi) * _q
    del _q, phi

    print('Qx: shape = {}\n{}'.format(qx.shape, qx.reshape((25, 1))))

    return qx, qy, twotheta


# TODO FIXME - Review and remove
def next_test_q_resolution_weiren(gpsans_f):

    filename = gpsans_f['sample_scattering_2']
    ws = LoadHFIRSANS(Filename=filename)

    # logs
    sl = SampleLogs(ws)
    wavelength = sl.single_value("wavelength")
    wavelength_spread = sl.single_value("wavelength-spread")
    r1 = sl.single_value("sample-aperture-diameter") * 1e-3 / 2  # radius in m
    r2 = sl.single_value("source-aperture-diameter") * 1e-3 / 2  # radius in m
    # this is hard-coded in the resolution (bar scan no implementend yet)
    pixel_size_x = 0.0055  # m
    pixel_size_y = 0.0043  # m

    # Those are ordered by workspace index.
    qx_arr, qy_arr, dqx_arr, dqy_arr = calculate_q_dq(ws)

    spec_info = ws.spectrumInfo()
    l1 = spec_info.l1()
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            two_theta = spec_info.twoTheta(i)
            l2 = spec_info.l2(i)
            qx = qx_arr[i]
            qy = qy_arr[i]

            res_dqx, res_dqy = sigma_neutron_weiren(
                wavelength, wavelength_spread,
                qx, qy, two_theta,
                l1, l2,
                r1, r2,
                pixel_size_x, pixel_size_y,
            )

            print('Spectrum {} of {}  Ricardo = {}, Weiren = {}'
                  ''.format(i, spec_info.size(), dqx_arr[i], res_dqx))

            assert dqx_arr[i] == res_dqx
            assert dqy_arr[i] == res_dqy
