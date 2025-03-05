import numpy as np
import pytest

from drtsans.dataobjects import IQazimuthal
from drtsans.iq import BinningMethod, BinningParams, bin_annular_into_q1d
from tests.unit.drtsans.i_of_q_binning_tests_data import generate_test_data, generate_test_data_wavelength


# This module supports testing data for issue #246.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/246

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data


def test_1d_annular_no_wt():
    """Test annular binning I(Qx, Qy) with no-weight binning method

    Returns
    -------
    None
    """
    # Initialize range of phi angle and Q
    phi_min = 0
    phi_max = 360.0
    num_bins = 10

    q_min = 0.003
    q_max = 0.006

    # Generate testing data: Get Q2D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test the high level method
    # Define input data
    test_i_q = IQazimuthal(
        intensity=intensities,
        error=sigmas,
        qx=qx_array,
        qy=qy_array,
        delta_qx=dqx_array,
        delta_qy=dqy_array,
    )

    # Annular binning
    phi_binning = BinningParams(phi_min, phi_max, num_bins)
    binned_iq = bin_annular_into_q1d(test_i_q, phi_binning, q_min, q_max, BinningMethod.NOWEIGHT)

    # Verify I(Q), sigma I(Q) and dQ
    assert binned_iq.intensity[1] == pytest.approx(63.66666667, abs=1e-8), "Binned intensity is wrong"
    assert binned_iq.error[1] == pytest.approx(3.257470048, abs=1e-8), "Binned sigma I is wrong"
    np.testing.assert_almost_equal(binned_iq.phi, np.linspace(start=18, stop=phi_max - 18, num=num_bins))


def test_1d_annular_out_of_range_angles():
    """Test annular binning I(Qx, Qy) supplying the azimuthal angle outside of 0<azimuthal<360deg

    Returns
    -------
    None
    """
    # Initialize range of phi angle and Q
    phi_min = -90
    phi_max = 270.0
    num_bins = 10

    q_min = 0.003
    q_max = 0.006

    # Generate testing data: Get Q2D data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)

    # Test the high level method
    # Define input data
    test_i_q = IQazimuthal(
        intensity=intensities,
        error=sigmas,
        qx=qx_array,
        qy=qy_array,
        delta_qx=dqx_array,
        delta_qy=dqy_array,
    )

    # Annular binning
    phi_binning = BinningParams(phi_min, phi_max, num_bins)
    with pytest.raises(ValueError):
        bin_annular_into_q1d(test_i_q, phi_binning, q_min, q_max, BinningMethod.NOWEIGHT)


def test_1d_flat_data():
    # 2d array of ones with -10 <= Qx/Qy <= 10
    data2d = IQazimuthal(
        intensity=np.full((21, 21), 1.0, dtype=float),
        error=np.full((21, 21), 1.0, dtype=float),
        qx=np.arange(-10, 11, 1, dtype=float),
        qy=np.arange(-10, 11, 1, dtype=float),
    )

    # perform the annular binning
    phi_binning = BinningParams(0, 360, 18)  # 20 deg bins
    binned_iq = bin_annular_into_q1d(data2d, phi_binning, q_min=9.0, q_max=10.0, method=BinningMethod.NOWEIGHT)

    # verify the results
    assert binned_iq.intensity.size == 18
    np.testing.assert_equal(binned_iq.intensity, 1.0)
    # uncertainties are proportional to the number of input bins contributing to the annular wedge
    np.testing.assert_allclose(binned_iq.error, 0.6, atol=0.11)
    # actually the azimuthal angle
    np.testing.assert_equal(binned_iq.phi, np.arange(10.0, 360.0, 20.0))


def test_1d_flat_data_wl():
    """Test with and without wavelength binning with constant intensity"""
    # Q2D data for one wavelength
    data2d = IQazimuthal(
        intensity=np.full((21, 21), 1.0, dtype=float),
        error=np.full((21, 21), 1.0, dtype=float),
        qx=np.arange(-10, 11, 1, dtype=float),
        qy=np.arange(-10, 11, 1, dtype=float),
    )
    data2d = data2d.ravel()

    # repeat Q2D data 3 times for 3 wavelengths
    wl_uniq = [1.5, 2.5, 3.5]
    num_wl = len(wl_uniq)
    data2d_wl = IQazimuthal(
        intensity=np.tile(data2d.intensity, num_wl),
        error=np.tile(data2d.error, num_wl),
        qx=np.tile(data2d.qx, num_wl),
        qy=np.tile(data2d.qy, num_wl),
        wavelength=np.repeat(wl_uniq, len(data2d.intensity)),  # [1.5, 1.5, ..., 2.5, 2.5, ..., 3.5, 3.5, ...]
    )

    # perform the annular binning
    phi_binning = BinningParams(0, 360, 18)  # 20 deg bins
    binned_i1d_wl = bin_annular_into_q1d(
        data2d_wl, phi_binning, q_min=9.0, q_max=10.0, method=BinningMethod.NOWEIGHT, wavelength_bins=None
    )
    binned_i1d_no_wl = bin_annular_into_q1d(
        data2d_wl, phi_binning, q_min=9.0, q_max=10.0, method=BinningMethod.NOWEIGHT
    )

    # verify the results
    assert binned_i1d_wl.intensity.shape[0] == phi_binning.bins * num_wl
    assert binned_i1d_no_wl.intensity.shape[0] == phi_binning.bins
    # intensity
    np.testing.assert_equal(binned_i1d_wl.intensity, 1.0)
    np.testing.assert_equal(binned_i1d_no_wl.intensity, 1.0)
    # uncertainties are inversely proportional to the number of input bins contributing to the annular wedge
    np.testing.assert_allclose(binned_i1d_wl.error, 0.6, atol=0.11)
    np.testing.assert_allclose(binned_i1d_no_wl.error, 0.4, atol=0.12)
    # actually the azimuthal angle
    np.testing.assert_equal(binned_i1d_wl.mod_q, np.tile(np.arange(10.0, 360.0, 20.0), num_wl))


def test_1d_wavelengths():
    """Test no wavelength binning with more realistic intensity data"""
    # Define input Q2D data
    num_wl = 3
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array, wl_array = generate_test_data_wavelength(2, num_wl)
    data2d = IQazimuthal(
        intensity=intensities,
        error=sigmas,
        qx=qx_array,
        qy=qy_array,
        delta_qx=dqx_array,
        delta_qy=dqy_array,
        wavelength=wl_array,
    )

    # perform the annular binning
    phi_binning = BinningParams(0, 360, 18)  # 20 deg bins
    binned_i1d_wl = bin_annular_into_q1d(data2d, phi_binning, method=BinningMethod.NOWEIGHT, wavelength_bins=None)

    # verify the results
    assert binned_i1d_wl.intensity.shape[0] == phi_binning.bins * num_wl
    # intensity
    assert np.nanmin(binned_i1d_wl.intensity) == pytest.approx(39.333, rel=1e-4)
    assert np.nanmax(binned_i1d_wl.intensity) == pytest.approx(77.0, rel=1e-4)
    # error
    assert np.nanmin(binned_i1d_wl.error) == pytest.approx(2.9155, rel=1e-4)
    assert np.nanmax(binned_i1d_wl.error) == pytest.approx(5.0662, rel=1e-4)
    # actually the azimuthal angle
    np.testing.assert_equal(binned_i1d_wl.mod_q, np.tile(np.arange(10.0, 360.0, 20.0), num_wl))

    # verify one bin
    assert binned_i1d_wl.mod_q[3] == pytest.approx(70.0, rel=1e-4)
    assert binned_i1d_wl.intensity[3] == pytest.approx(59.0, rel=1e-4)
    assert binned_i1d_wl.error[3] == pytest.approx(4.4347, rel=1e-4)
    assert binned_i1d_wl.wavelength[3] == pytest.approx(1.5, rel=1e-4)


if __name__ == "__main__":
    pytest.main([__file__])
