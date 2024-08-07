from drtsans.plots import plot_IQmod, plot_IQazimuthal, plot_detector
from drtsans.plots.api import _save_file, MatBackendManager, Backend
from drtsans.dataobjects import IQmod, IQazimuthal
from mantid.simpleapi import LoadEmptyInstrument, LoadNexus
import numpy as np
import os
import pytest
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread
from typing import Tuple, Any
import mpld3  # noqa E402


def verify_images(test_png: str, gold_png):
    """Verify the image output from test is same as what is expected (gold)

    AssertionError will be raised if they do not match
    """
    # check
    assert os.path.exists(gold_png), f"Gold/reference data file {gold_png} cannot be located"
    assert os.path.exists(test_png), f"Test result data file {test_png} cannot be located"

    # Import PNG to 2d array
    tested_image = imread(test_png)
    gold_image = imread(gold_png)

    # Verify
    np.testing.assert_allclose(
        tested_image,
        gold_image,
        err_msg=f"Testing result {tested_image} does not match " f"the expected result {gold_image}",
    )


def fileCheckAndRemove(filename, remove=True):
    """Convenience function for doing simple checks that the file was created.
    The ``remove`` option is available to make debugging new tests easier."""
    assert os.path.exists(filename), '"{}" does not exist'.format(filename)
    assert os.path.getsize(filename) > 100, '"{}" is too small'.format(filename)
    if remove:
        os.remove(filename)


@pytest.mark.parametrize(
    "backend, filename",
    [("mpl", "test_IQmod.png"), ("d3", "test_IQmod.json")],
    ids=["mpl", "d3"],
)
def test_IQmod(backend, filename):
    """Test plotting single a IQmod dataset"""
    x = np.linspace(0.0, 4 * np.pi, 50)
    e = np.zeros(50) + 0.1
    data = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    plot_IQmod([data], filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename)


@pytest.mark.parametrize(
    "backend, filename",
    [("mpl", "test_IQmod_multi.png"), ("d3", "test_IQmod_multi.json")],
    ids=["mpl", "d3"],
)
def test_IQmod_multi(backend, filename):
    """Test over-plotting multiple IQmod datasets"""
    x = np.linspace(0.0, 4 * np.pi, 50)
    e = np.zeros(50) + 0.1
    data1 = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    data2 = IQmod(intensity=np.cos(x), error=e, mod_q=x)

    plot_IQmod([data1, data2], filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename)


@pytest.mark.parametrize(
    "backend, filename, show",
    [
        ("mpl", "test_IQmod_multi.png", False),
        ("mpl", "test_IQmod_multi.png", True),
        ("d3", "test_IQmod_multi.json", False),
        ("mpld3", "test_IQmod_multi.json", False),
    ],
    ids=["mpl", "mpl", "d3", "mpld3"],
)
def test_save_file(backend, filename, show):
    """Test save_file"""

    figure, _ax = plt.subplots()

    def mock_show():
        return False

    figure.show = mock_show

    _save_file(figure, filename, Backend.getMode(backend), show)
    fileCheckAndRemove(filename)


def test_save_file_inline_error():
    """Test save_file with unsupported backend"""

    figure, _ax = plt.subplots()

    def mock_show(figure):
        return False

    mpld3.show = mock_show
    backend, filename, show = "inline", "test_IQmod_multi.json", True
    with pytest.raises(RuntimeError) as excinfo:
        _save_file(figure, filename, Backend.getMode(backend), show)
    assert str(excinfo.value) == "Unsupported backend: inline"


def test_MatBackendManager():
    """Test MatBackendManager"""

    outer_backend = "pdf"
    inner_backend = "agg"

    matplotlib.use(outer_backend)
    with MatBackendManager(inner_backend) as _m:
        assert matplotlib.get_backend() == inner_backend
    assert matplotlib.get_backend() == outer_backend


@pytest.fixture
def test_iq2d_data() -> Tuple[Any, Any, Any, Any]:
    """Generate Qx, Qy, I(qx, qy), Error(qx, qy) data"""
    # Qx: 60 values, Qy: 40 values
    x = np.linspace(0.0, 4 * np.pi, 60) - 3.0
    y = np.linspace(0.5 * np.pi, 4.5 * np.pi, 40) - 3.0
    # Calculate intensity and error
    mesh_x, mesh_y = np.meshgrid(x, y, sparse=False, indexing="ij")
    intensity = np.abs(np.sin(mesh_x) + np.cos(mesh_y))
    error = np.sqrt(intensity)

    # Transfer to correct orientation
    mesh_x = mesh_x.T
    mesh_y = mesh_y.T
    intensity = intensity.T
    error = error.T

    # now try to find zero...
    x_zero_index = np.argmin(np.abs(x))
    y_zero_index = np.argmin(np.abs(y))
    # mask center
    for x_index in [-1, 0, 1]:
        for y_index in [-1, 0, 1]:
            center_x_index = x_zero_index + x_index
            center_y_index = y_zero_index + y_index
            intensity[center_y_index][center_x_index] = np.nan
            error[center_y_index][center_x_index] = np.nan

    # mark lower right corner
    intensity[1, 59] = 8.0

    assert intensity.shape == (40, 60), (
        f"Expected intensity is 40 row (Qy) and 60 column (Qx) " f"but not {intensity.shape}"
    )

    return mesh_x, mesh_y, intensity, error


@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_1d.png", "tests/unit/references/gold_IQazimuthal_2d_T.png"),
        ("d3", "test_IQazimuthal_1d.json", None),
    ],
    ids=["mpl_test", "d3"],
)
def test_IQazimuthal_1d(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 1d Qx and Qy arrays"""
    # Generate input arrays
    mesh_x, mesh_y, intensity, error = test_iq2d_data
    x = mesh_x[0]
    assert x.min() < x.max()
    y = mesh_y[:, 0]
    assert y.min() < y.max()
    # construct IQazimuthal: following bin_iq_2d() routine, i.e., intensity is (num_qx, num_qy)
    assert intensity.T.shape[0] == len(x)
    assert intensity.T.shape[1] == len(y)
    data = IQazimuthal(intensity=intensity.T, error=error.T, qx=x, qy=y)

    # plot
    plot_IQazimuthal(data, filename=filename, backend=backend)
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=False)


@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_2d.png", "tests/unit/references/gold_IQazimuthal_2d_T.png"),
        ("d3", "test_IQazimuthal_2d.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity, error=error, qx=x, qy=y)

    # plot
    plot_IQazimuthal(data, filename=filename, backend=backend)
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_2d_T.png", "tests/unit/references/gold_IQazimuthal_2d_T.png"),
        ("d3", "test_IQazimuthal_2d_T.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d_transposed(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity.T, error=error.T, qx=x.T, qy=y.T)

    # plot
    plot_IQazimuthal(data, filename=filename, backend=backend)
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


# Broken!
@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_2d_selections.png", "tests/unit/references/new_IQazimuthal_2d_selections.png"),
        ("d3", "test_IQazimuthal_2d_selections.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d_selections(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity, error=error, qx=x, qy=y)

    # plot
    plot_IQazimuthal(
        data,
        filename=filename,
        backend=backend,
        qmin=0.0,
        qmax=9.0,
        wedges=((-30.0, 30.0),),
    )
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


# Broken
@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        (
            "mpl",
            "test_IQazimuthal_2d_asymmetric_wedge.png",
            "tests/unit/references/new_IQazimuthal_2d_asymmetric_wedge.png",
        ),
        ("d3", "test_IQazimuthal_2d_asymmetric_wedge.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d_asymmetric_wedge(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity.T, error=error.T, qx=x.T, qy=y.T)

    # plot
    plot_IQazimuthal(
        data,
        filename=filename,
        backend=backend,
        wedges=[(-35.0, 45.0), (125.0, 145.0)],
        symmetric_wedges=False,
    )
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_2d_ring.png", "tests/unit/references/gold_IQazimuthal_2d_ring.png"),
        ("d3", "test_IQazimuthal_2d_ring.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d_ring(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity, error=error, qx=x, qy=y)

    # plot
    plot_IQazimuthal(
        data,
        filename=filename,
        backend=backend,
        qmin=1.0,
        qmax=2.0,
    )
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


@pytest.mark.parametrize(
    "backend, filename, reference_name",
    [
        ("mpl", "test_IQazimuthal_2d_ring_T.png", "tests/unit/references/gold_IQazimuthal_2d_ring.png"),
        ("d3", "test_IQazimuthal_2d_ring_T.json", None),
    ],
    ids=["mpl", "d3"],
)
def test_IQazimuthal_2d_transposed_ring(backend, filename, reference_name, test_iq2d_data):
    """Test plotting IQazimuthal with 2d Qx and Qy arrays"""
    # construct IQazimuthal
    x, y, intensity, error = test_iq2d_data
    data = IQazimuthal(intensity=intensity.T, error=error.T, qx=x.T, qy=y.T)

    # plot
    plot_IQazimuthal(data, filename=filename, backend=backend, qmin=1.0, qmax=2.0)
    plt.close()

    # verify
    if reference_name:
        fileCheckAndRemove(filename, remove=False)
        verify_images(filename, reference_name)
    fileCheckAndRemove(filename, remove=True)


@pytest.mark.parametrize(
    "backend, filename",
    [("mpl", "test_detector.png"), ("d3", "test_detector.json")],
    ids=["mpl", "d3"],
)
def test_detector(backend, filename):
    """Test plotting in detector space from a mantid workspace"""
    workspace = LoadEmptyInstrument(InstrumentName="CG3")  # this will load monitors as well
    plot_detector(workspace, filename, backend)
    plt.close()
    fileCheckAndRemove(filename, False)


@pytest.mark.datarepo
def test_xaxis_direction(datarepo_dir, clean_workspace):
    r"""Test values of X-axis in plot_detector decrease when looking at the picture from left to right"""
    # wing_detector.nxs contains intensities for the wing detector that can be plotted
    workspace = LoadNexus(os.path.join(datarepo_dir.sans, "plots", "wing_detector.nxs"))
    clean_workspace(workspace)
    filename = "test_xaxis_direction.png"
    plot_detector(
        workspace,
        filename=filename,
        backend="mpl",
        panel_name="wing_detector",
        axes_mode="xy",
    )
    plt.close()
    fileCheckAndRemove(filename, remove=True)


if __name__ == "__main__":
    pytest.main([__file__])
