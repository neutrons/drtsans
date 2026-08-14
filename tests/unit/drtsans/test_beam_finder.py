# local imports
from drtsans.beam_finder import BeamCenterNotFound, _calculate_neutron_drop, fbc_options_json, find_beam_center
from drtsans.mono.biosans.beam_finder import find_beam_center as biosans_find_beam_center
from drtsans.redparams import reduction_parameters

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
    assert mock_logger.warning.call_count == 1
    assert "not finite" in mock_logger.warning.call_args.args[0]
    # a substitution is a caveat, not a failure: an error would fail the whole autoreduction
    assert mock_logger.error.call_count == 0


def test_find_beam_center_fallback_only_for_one_axis(fitted_center):
    """A fallback for the failing axis is enough; the other axis keeps its fitted value"""
    mock_logger = fitted_center((float("nan"), 0.3))

    x, y, center_type, _ = find_beam_center("unused", fallback_center=(FALLBACK[0], None))

    assert (x, y) == pytest.approx((FALLBACK[0], 0.3))
    assert center_type == "fallback"
    assert mock_logger.warning.call_count == 1
    assert mock_logger.error.call_count == 0


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


@pytest.mark.parametrize(
    "instrument_name, expected",
    [("EQSANS", FALLBACK), ("GPSANS", (0.0, 0.0)), ("BIOSANS", (0.0, 0.0))],
)
def test_fbc_options_json_instrument_defaults(instrument_name, expected):
    """Opting in hands find_beam_center the coordinates carried by the instrument's own schema"""
    parameters = reduction_parameters(
        {
            "instrumentName": instrument_name,
            "iptsNumber": 42,
            "beamCenter": {"runNumber": 1, "useFallbackBeamCenter": True},
        },
        validate=False,
    )

    assert fbc_options_json(parameters)["fallback_center"] == pytest.approx(expected)


@pytest.mark.parametrize("fallback_center", [None, FALLBACK])
def test_biosans_find_beam_center_forwards_fallback(mocker, fallback_center):
    """BIOSANS derives the wing and midrange centers from the main detector, so it must relay the fallback"""
    mocker.patch("drtsans.mono.biosans.beam_finder.mtd")
    mocker.patch("drtsans.mono.biosans.beam_finder._beam_center_gravitational_drop", return_value=0.1)
    mock_find = mocker.patch(
        "drtsans.mono.biosans.beam_finder.bf.find_beam_center",
        return_value=(0.5, 0.3, "calculated", {}),
    )
    # non-zero distances keep the instrument geometry out of this test
    distances = dict(
        sample_det_cent_main_detector=1.0,
        sample_det_cent_wing_detector=1.1,
        sample_det_cent_midrange_detector=1.2,
    )
    options = {} if fallback_center is None else {"fallback_center": fallback_center}

    biosans_find_beam_center("unused", **distances, **options)

    # an omitted fallback must reach drtsans.find_beam_center as the fatal (None, None)
    expected = (None, None) if fallback_center is None else fallback_center
    assert mock_find.call_args.kwargs["fallback_center"] == expected


if __name__ == "__main__":
    pytest.main([__file__])
