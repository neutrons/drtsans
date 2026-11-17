# Unit tests for drtsans.auto_wedge
import numpy as np
import pytest

from drtsans.auto_wedge import SIGMA_TO_FWHM, PeakFit, RingStatus, _circular_difference


class TestRingSelection:
    """Building blocks that select the Q-rings used to determine the wedge angles"""

    @staticmethod
    def _peak(**changes) -> PeakFit:
        """Peak at 181 deg with height 0.2 +/- 0.002 and FWHM 17 +/- 0.5 deg, with optional field changes"""
        values = dict(height=0.2, height_error=0.002, centre=181.0, centre_error=0.1, fwhm=17.0, fwhm_error=0.5)
        values.update(changes)
        return PeakFit(**values)

    @pytest.mark.parametrize(
        "angle, reference, expected",
        [(361.0, 1.0, 0.0), (1.0, 359.0, 2.0), (359.0, 1.0, -2.0), (181.0, 1.0, -180.0), (90.0, 450.0, 0.0)],
    )
    def test_circular_difference(self, angle, reference, expected):
        assert _circular_difference(angle, reference) == pytest.approx(expected)

    def test_ring_status_accepted(self):
        accepted = {status for status in RingStatus if status.accepted}
        assert accepted == {RingStatus.SEED, RingStatus.RESCUED}

    def test_peak_fit_from_fit_result(self):
        # Arrange: a negative Sigma describes the same peak, since the Gaussian depends only on Sigma squared
        result = {"f2.Height": (0.2, 0.002), "f2.PeakCentre": (361.3, 0.1), "f2.Sigma": (-7.5, 0.2)}
        # Act
        peak = PeakFit.from_fit_result(result, "f2")
        # Assert
        assert (peak.height, peak.height_error, peak.centre, peak.centre_error) == (0.2, 0.002, 361.3, 0.1)
        assert peak.fwhm == pytest.approx(7.5 * SIGMA_TO_FWHM)
        assert peak.fwhm_error == pytest.approx(0.2 * SIGMA_TO_FWHM)

    def test_peak_fit_significance_and_contrast(self):
        peak = self._peak()
        assert peak.significance == pytest.approx(100.0)
        assert peak.contrast(0.25) == pytest.approx(0.8)
        for background in (0.0, -0.1, np.nan, np.inf):
            assert np.isnan(peak.contrast(background))
        for height_error in (0.0, -0.002, np.nan):
            assert np.isnan(self._peak(height_error=height_error).significance)

    @pytest.mark.parametrize(
        "changes", [{}, {"height_error": 0.0}, {"centre_error": -0.1}, {"fwhm_error": np.nan}, {"centre": np.inf}]
    )
    def test_peak_fit_is_valid(self, changes):
        assert self._peak(**changes).is_valid == (not changes)

    def test_peak_fit_weights_and_with_centre(self):
        peak = self._peak()
        moved = peak.with_centre(1.3)
        assert peak.centre_weight == pytest.approx(10.0)
        assert peak.fwhm_weight == pytest.approx(2.0)
        assert moved.centre == 1.3 and peak.centre == 181.0  # the original is unchanged
        assert moved.fwhm == peak.fwhm and moved.height == peak.height
