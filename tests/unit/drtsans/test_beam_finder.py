# local imports
from drtsans.beam_finder import BeamCenterNotFound, _calculate_neutron_drop, find_beam_center

# third party imports
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# standard imports


# EQSANS beam center assumed when a fit fails, in meters
FALLBACK = (0.025239, 0.0170801)

NON_FINITE = [
    (float("nan"), float("nan")),
    (float("nan"), 0.3),
    (0.5, float("nan")),
    (float("inf"), 0.3),
]


def test_calculate_neutron_drop():
    path_length = 15.0  # meters
    wavelength = 18.0  # Angstrom
    gravity_drop = _calculate_neutron_drop(path_length, wavelength)
    assert_almost_equal(gravity_drop, 0.02284, decimal=4)


@pytest.fixture
def fitted_center(mocker):
    """Make find_beam_center return a chosen (X, Y) without needing a workspace"""

    def _fitted_center(coordinates):
        mocker.patch("drtsans.beam_finder.Integration")
        mocker.patch("drtsans.beam_finder.mask_spectra_with_special_values")
        mocker.patch("drtsans.beam_finder.solid_angle_correction")
        mocker.patch("drtsans.beam_finder.DeleteWorkspace")
        mocker.patch("drtsans.beam_finder.FindCenterOfMassPosition", return_value=coordinates)
        return mocker.patch("drtsans.beam_finder.logger")

    return _fitted_center


@pytest.mark.parametrize("coordinates", NON_FINITE)
def test_find_beam_center_no_fallback_raises(fitted_center, coordinates):
    """Without a fallback, a non-finite coordinate stops the reduction at its source"""
    mock_logger = fitted_center(coordinates)

    with pytest.raises(BeamCenterNotFound, match="not finite"):
        find_beam_center("unused")

    # the exception message carries the report, so it is not logged a second time
    assert mock_logger.error.call_count == 0


@pytest.mark.parametrize("coordinates", NON_FINITE)
def test_find_beam_center_fallback_substituted(fitted_center, coordinates):
    """With a fallback, the failing coordinate is assumed and reported in one message"""
    mock_logger = fitted_center(coordinates)

    x, y, center_type, _ = find_beam_center("unused", fallback_center=FALLBACK)

    for axis, value, fitted, fallback in zip("XY", (x, y), coordinates, FALLBACK):
        expected = fitted if np.isfinite(fitted) else fallback
        assert value == pytest.approx(expected), f"unexpected {axis}"
    assert center_type == "fallback"
    # a single composite message, not one per offending axis
    assert mock_logger.error.call_count == 1
    assert "not finite" in mock_logger.error.call_args.args[0]


def test_find_beam_center_fallback_only_for_one_axis(fitted_center):
    """A fallback for the failing axis is enough; the other axis keeps its fitted value"""
    mock_logger = fitted_center((float("nan"), 0.3))

    x, y, center_type, _ = find_beam_center("unused", fallback_center=(FALLBACK[0], None))

    assert (x, y) == pytest.approx((FALLBACK[0], 0.3))
    assert center_type == "fallback"
    assert mock_logger.error.call_count == 1


def test_find_beam_center_fallback_missing_for_failing_axis(fitted_center):
    """A fallback for the healthy axis does not rescue the failing one"""
    fitted_center((0.5, float("nan")))

    with pytest.raises(BeamCenterNotFound, match="no fallback value is available for Y"):
        find_beam_center("unused", fallback_center=(FALLBACK[0], None))


@pytest.mark.parametrize("fallback_center", [(None, None), FALLBACK])
def test_find_beam_center_finite_untouched(fitted_center, fallback_center):
    """A converged fit passes through whatever the fallback is"""
    mock_logger = fitted_center((0.5, 0.3))

    x, y, center_type, _ = find_beam_center("unused", fallback_center=fallback_center)

    assert (x, y) == pytest.approx((0.5, 0.3))
    assert center_type == "calculated"
    assert mock_logger.error.call_count == 0


if __name__ == "__main__":
    pytest.main([__file__])
