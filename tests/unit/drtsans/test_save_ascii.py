from drtsans.dataobjects import IQazimuthal, I1DAnnular
from drtsans.save_ascii import load_ascii_binned_2D, save_ascii_binned_2D, save_ascii_binned_1D_annular
import numpy as np
from pathlib import Path
import pytest
from tempfile import NamedTemporaryFile


def createIQAzimuthal(with_resolution=True):
    """Create :pyobj:`~drtsans.dataobjects.IQazimuthal`. The data has shape=(21,21)

    Parameters
    ----------
    with_resolution: bool
        Whether or not to add the resolution fields
    """
    resolutionArgs = {}
    if with_resolution:
        resolutionArgs["delta_qx"] = np.full((21, 21), 1.0, dtype=float)
        resolutionArgs["delta_qy"] = np.full((21, 21), 1.0, dtype=float)

    return IQazimuthal(
        intensity=np.arange(21 * 21, dtype=float).reshape((21, 21)),
        error=np.full((21, 21), 1.0, dtype=float),
        qx=np.arange(-10, 11, 1, dtype=float),
        qy=np.arange(-10, 11, 1, dtype=float),
        **resolutionArgs,
    )


def assert_IQazi_allclose(data2d_exp, data2d_obs, err_msg=""):
    if err_msg:
        err_msg += " "
    for key in ["intensity", "error", "qx", "qy", "delta_qx", "delta_qy", "wavelength"]:
        exp = getattr(data2d_exp, key)
        obs = getattr(data2d_obs, key)
        if exp is None or obs is None:
            assert exp == obs, "{} are not both None".format(key)
        else:
            np.testing.assert_allclose(exp, obs, err_msg="{}{} doesn't match".format(err_msg, key))


@pytest.mark.parametrize("with_resolution", [True, False])
def test_ascii_binned_2D_roundtrip(with_resolution):
    # create test data with -10 <= Qx/Qy <= 10
    data2d = createIQAzimuthal(with_resolution)

    # file to write to
    filename = Path(NamedTemporaryFile(suffix=".dat").name)

    # write out the data
    print("writing data to", filename)
    save_ascii_binned_2D(filename, "test data", data2d)
    assert filename.exists()

    # read the data back in
    print("reading data from", filename)
    data2d_reread = load_ascii_binned_2D(filename)
    assert data2d_reread

    # flatten the data objects for comparison
    data2d = data2d.ravel()
    data2d_reread = data2d_reread.ravel()

    # validate the data is unchanged
    assert_IQazi_allclose(data2d, data2d_reread)

    # cleanup - not using fixture so failures can be inspected
    filename.unlink()


def assert_I1DAnnular_allclose(data1d_exp, data1d_obs, err_msg=""):
    if err_msg:
        err_msg += " "
    for key in ["intensity", "error", "phi"]:
        exp = getattr(data1d_exp, key)
        obs = getattr(data1d_obs, key)
        if exp is None or obs is None:
            assert exp == obs, "{} are not both None".format(key)
        else:
            np.testing.assert_allclose(exp, obs, err_msg="{}{} doesn't match".format(err_msg, key))


def test_save_ascii_annular():
    # create test data
    data1d = I1DAnnular(
        intensity=[1, 2, 3],
        error=[4, 5, 6],
        phi=[7, 8, 9],
    )

    # file to write to
    filename = Path(NamedTemporaryFile(suffix=".dat").name)

    # write out the data
    print("writing data to", filename)
    save_ascii_binned_1D_annular(filename, "test data", data1d)
    assert filename.exists()

    # read the data back in
    print("reading data from", filename)
    csv_data = np.genfromtxt(filename, comments="#", dtype=np.float64, skip_header=2)
    num_cols = len(csv_data[0])
    assert num_cols == 3, "Incompatible number of columns: should be 3"
    data1d_reread = I1DAnnular(
        phi=csv_data[:, 0],
        intensity=csv_data[:, 1],
        error=csv_data[:, 2],
    )  # wavelength isn't in the file
    del csv_data

    # validate the data is unchanged
    assert_I1DAnnular_allclose(data1d, data1d_reread)


if __name__ == "__main__":
    pytest.main([__file__])
