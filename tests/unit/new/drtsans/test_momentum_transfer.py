import pytest
import numpy as np
from pytest import approx
from drtsans.geometry import pixel_centers
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.momentum_transfer import convert_to_q, _filter_and_replicate, subpixel_info
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


@pytest.fixture(scope='module')
@namedtuplefy
def data_subpixel_info():
    r"""Data to be used for `test_subpixel_info`"""
    return dict(number_pixels=4,
                wavelength_bin_boundaries=[2.0, 2.1],
                intensities=[[1.0, 2.0],  # first tube
                             [3.0, 2.0]],  # second tube
                # Pixel (x, y) coordinates, in mili-meters
                pixel_positions=[[10, -10], [10, 10], [-10, -10], [-10, 10]],
                n_horizontal=2,  # number of subpixel along the X-axis
                n_vertical=4,  # number of subpixel along the Y-axis
                number_subpixels=8,  # just n_horizontal * n_vertical
                # subpixel coordinates for a pixel centered at (x, y) == (0, 0)
                subpixel_positions=[[5, -7.5], [5, -2.5], [5, 2.5], [5, 7.5],
                                    [-5, -7.5], [-5, -2.5], [-5, 2.5], [-5, 7.5]],
                sample_detector_distance=1.0  # all subpixels have same Z coordinate
                )


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack',
                           'n_tubes': 2, 'n_pixels': 2, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 1.0,  # detector center
                           'diameter': 0.02, 'height': 0.02}],  # pixel dimensions
                         indirect=True)
def test_subpixel_info(data_subpixel_info, workspace_with_instrument):
    r"""
    Subpixel binning for a detector array made up of two tubes, each with two pixels.
    There's no spacing between tubes.
    Detector is 1 meter away from the sample
    The shape of a detector pixel is a cylinder of 20mm diameter and 20mm in height.
    """
    data = data_subpixel_info  # handy shortcut
    input_workspace = unique_workspace_dundername()  # temporary workspace
    workspace = workspace_with_instrument(axis_values=data.wavelength_bin_boundaries,
                                          intensities=data.intensities,
                                          output_workspace=input_workspace)
    workspace.detectorInfo().setMasked(2, True)  # Mask the third pixel
    spectrum_info = workspace.spectrumInfo()
    get_spectrum_definition = spectrum_info.getSpectrumDefinition
    info_indexes = [get_spectrum_definition(idx)[0][0] for idx in range(workspace.getNumberHistograms())]
    # test the positions of the pixels against expected_pixel_positions
    pixel_positions = pixel_centers(input_workspace, info_indexes)
    assert 1.e03 * pixel_positions[:, :-1] == pytest.approx(np.array(data.pixel_positions))
    # Find the position of the subpixels in polar coordinates (info data structure)
    info = subpixel_info(input_workspace, data.n_horizontal, data.n_vertical)
    # Transform from polar to cartesian coordinates
    x = info.l2 * np.sin(info.two_theta) * np.cos(info.azimuthal)
    y = info.l2 * np.sin(info.two_theta) * np.sin(info.azimuthal)
    z = info.l2 * np.cos(info.two_theta)
    assert z == pytest.approx(data.sample_detector_distance * np.ones(len(info.keep) * data.number_subpixels))
    # array `xy` contains the coordinates of the subpixels, in mili-meters
    # xy.shape = (number_pixels, number_subpixels, 2)
    xy = np.column_stack((x, y)).reshape((len(info.keep), data.number_subpixels, 2))
    # Translate the subpixels by removing the center pixel positions, so that each center pixel is at (0, 0)
    # After translation, all four set of subpixels have the same coordinates
    xy = 1000 * (xy - pixel_positions[info.keep, :-1][:, None, :])
    for subpixel_positions in xy:  # unmasked pixels
        assert subpixel_positions == pytest.approx(np.array(data.subpixel_positions))


def test_filter_and_replicate():
    def to_np(alist):
        return np.array(alist).reshape((len(alist), 1))
    lam, intensity, error = to_np([1., 2., 3.]), to_np([100., 81., 64.]), to_np([10., 9., 8.])
    keep = np.array([1, 2])
    lam, intensity, error = _filter_and_replicate((lam, intensity, error), keep, 2, 1)
    assert lam == pytest.approx([2., 2., 3., 3.])
    assert intensity == pytest.approx([81., 81., 64., 64.])
    assert error == pytest.approx([9., 9., 8., 8.])


if __name__ == '__main__':
    pytest.main([__file__])
