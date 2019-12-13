import pytest
import numpy as np
from pytest import approx
from drtsans.momentum_transfer import convert_to_q, convert_to_subpixels
from drtsans.dataobjects import IQazimuthal
from mantid.simpleapi import MaskDetectors


def fake_resolution1(*args, **kwargs):
    """
    resolution dependent only on q
    """
    mode = kwargs.get('mode')
    if mode == 'scalar':
        q = args[0]
        return np.arange(len(q)).reshape(-1, 1)
    elif mode == 'azimuthal':
        qx = args[0]
        qy = args[1]
        return qx+qy, qx-qy
    else:
        raise NotImplementedError('fake_resolution1 mode')


def fake_resolution2(*args, **kwargs):
    """
    resolution dependent on q and geometry
    """
    mode = kwargs.get('mode')
    info = kwargs.get('pixel_info', None)
    if mode == 'scalar':
        q = args[0]
        if info:
            azi = info.azimuthal
            num_bins = q.shape[1]
            return np.repeat(azi, num_bins).reshape(-1, num_bins)
    else:
        raise NotImplementedError('fake_resolution1 mode')


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[10, 20],
                                           [30, 40]],
                           'uncertainties':[[1, 2],
                                            [3, 4]],
                           'axis_values': [5.9, 6.1],
                           'name':'BIOSANS', 'dx': 1, 'dy': 1}],
                         indirect=True)
def test_convert_to_mod_q(generic_workspace):
    ws = generic_workspace
    MaskDetectors(ws, WorkspaceIndexList=[2])
    intensity, error, modq, dq, lam = convert_to_q(ws, mode='scalar')
    assert intensity == approx([10, 20, 40], abs=1e-5)
    assert error == approx([1, 2, 4], abs=1e-5)
    assert dq == approx([0, 0, 0], abs=1e-5)
    assert lam == approx([6, 6, 6], abs=1e-5)
    # position of the detectors are at (0,5, 0.5, 5.0) and symmetric
    two_theta = np.arccos(5./np.sqrt(25.5))
    assert ws.spectrumInfo().twoTheta(0) == approx(two_theta, abs=1e-5)

    q = 4. * np.pi * np.sin(two_theta * 0.5) / 6.
    assert modq == approx([q, q, q], abs=1e-5)

    # test namedtuple result and resolution
    result = convert_to_q(ws, mode='scalar', resolution_function=fake_resolution1)
    assert result.delta_mod_q == approx([0, 1, 3], abs=1e-5)

    # test geometry dependent resolution
    result = convert_to_q(ws, mode='scalar', resolution_function=fake_resolution2)
    # note that x axis is pointing to the left, so the azimuthal angles by
    # spectrum number are -45, 45, -135, 135
    assert np.degrees(result.delta_mod_q) == approx([-45, 45, 135], abs=1e-5)


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[10, 20],
                                           [30, 40]],
                           'uncertainties':[[1, 2],
                                            [3, 4]],
                           'axis_values': [5.9, 6.1],
                           'name':'BIOSANS', 'dx': 1, 'dy': 1}],
                         indirect=True)
def test_convert_q_azimuthal(generic_workspace):
    ws = generic_workspace
    MaskDetectors(ws, WorkspaceIndexList=[2])
    result = convert_to_q(ws, mode='azimuthal')
    assert result.intensity == approx([10, 20, 40], abs=1e-5)
    assert result.error == approx([1, 2, 4], abs=1e-5)
    assert result.delta_qx == approx([0, 0, 0], abs=1e-5)
    assert result.delta_qy == approx([0, 0, 0], abs=1e-5)
    assert result.wavelength == approx([6, 6, 6], abs=1e-5)
    two_theta = np.arccos(5./np.sqrt(25.5))
    q = 4. * np.pi * np.sin(two_theta * 0.5) / 6.
    assert result.qx * (-1) == approx([q * np.sqrt(0.5), q*np.sqrt(0.5), -q*np.sqrt(0.5)], abs=1e-5)
    assert result.qy == approx([-q * np.sqrt(0.5), q*np.sqrt(0.5), q*np.sqrt(0.5)], abs=1e-5)

    # test with simple resolution
    result = convert_to_q(ws, mode='azimuthal', resolution_function=fake_resolution1)
    assert result.delta_qx == approx([-q * np.sqrt(2.), 0, q * np.sqrt(2.)], abs=1e-5)
    assert result.delta_qy == approx([0, -q * np.sqrt(2.), 0], abs=1e-5)


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[10, 20],
                                           [30, 40]],
                           'uncertainties':[[1, 2],
                                            [3, 4]],
                           'axis_values': [5.9, 6.1],
                           'name':'BIOSANS', 'dx': 1, 'dy': 1}],
                         indirect=True)
def test_convert_q_crystal(generic_workspace):
    ws = generic_workspace
    MaskDetectors(ws, WorkspaceIndexList=[2])
    result = convert_to_q(ws, mode='crystallographic')
    assert result.intensity == approx([10, 20, 40], abs=1e-5)
    assert result.error == approx([1, 2, 4], abs=1e-5)
    assert result.delta_qx == approx([0, 0, 0], abs=1e-5)
    assert result.delta_qy == approx([0, 0, 0], abs=1e-5)
    assert result.delta_qz == approx([0, 0, 0], abs=1e-5)
    assert result.wavelength == approx([6, 6, 6], abs=1e-5)
    two_theta = np.arccos(5./np.sqrt(25.5))
    kf = 2. * np.pi / 6.
    assert result.qx == approx([kf * np.sin(two_theta) * np.sqrt(0.5),
                                kf * np.sin(two_theta) * np.sqrt(0.5),
                                -kf * np.sin(two_theta) * np.sqrt(0.5)], abs=1e-5)
    assert result.qy == approx([-kf * np.sin(two_theta) * np.sqrt(0.5),
                                kf * np.sin(two_theta) * np.sqrt(0.5),
                                kf * np.sin(two_theta) * np.sqrt(0.5)], abs=1e-5)
    assert result.qz == approx([kf * (np.cos(two_theta) - 1),
                                kf * (np.cos(two_theta) - 1),
                                kf * (np.cos(two_theta) - 1)], abs=1e-5)

    # test with simple resolution - it should thow not implemented
    with pytest.raises(NotImplementedError):
        assert convert_to_q(ws, mode='crystallographic', resolution_function=fake_resolution1)


@pytest.mark.parametrize('generic_workspace',
                         [{'Nx': 2, 'Ny': 2,
                           'name': 'GPSANS', 'dx': 1, 'dy': 1}],
                         indirect=True)
def test_subpixel_binning(generic_workspace):
    ws = generic_workspace
    # create data
    qx = np.array([-0.001, -0.001, 0.001, 0.001])
    qy = np.array([-0.001, 0.001, -0.001, 0.001])
    intensity = np.array([52., 34., 90., 75.])
    sigma_Qx = np.array([15e-5, 13e-5, 16e-5, 14e-5])
    sigma_Qy = np.array([15e-5, 13e-5, 16e-5, 14e-5])
    i_q_azi = IQazimuthal(intensity, np.sqrt(intensity), qx, qy, delta_qx=sigma_Qx, delta_qy=sigma_Qy)
    # calculate subpixels
    res = convert_to_subpixels(i_q_azi, ws, 'azimuthal')
    # Q offsets from the original pixel
    d_qx = 0.0005
    d_qy = 0.0005
    # loop over the original pixels
    for i in range(4):
        assert res.intensity[i:16:4] == approx([intensity[i]] * 4)
        assert res.error[i:16:4] == approx([np.sqrt(intensity[i])] * 4)
        assert res.delta_qx[i:16:4] == approx([sigma_Qx[i]] * 4)
        assert res.delta_qy[i:16:4] == approx([sigma_Qy[i]] * 4)
        assert res.qx[i:16:4] == approx(qx[i] + np.array([-d_qx, -d_qx, d_qx, d_qx]))
        assert res.qy[i:16:4] == approx(qy[i] + np.array([-d_qy, d_qy, -d_qy, d_qy]))


if __name__ == '__main__':
    pytest.main([__file__])
