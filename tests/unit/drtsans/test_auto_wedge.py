# Unit tests for drtsans.auto_wedge
import numpy as np
import pytest

from drtsans.auto_wedge import SIGMA_TO_FWHM, PeakFit, RingFit, RingStatus, _circular_difference


@pytest.mark.parametrize(
    "angle, reference, expected",
    [(361.0, 1.0, 0.0), (1.0, 359.0, 2.0), (359.0, 1.0, -2.0), (181.0, 1.0, -180.0), (90.0, 450.0, 0.0)],
)
def test_circular_difference(angle, reference, expected):
    """Signed angle differences are wrapped into [-180, 180) deg"""
    assert _circular_difference(angle, reference) == pytest.approx(expected)


class TestRingStatus:
    """Outcome of the ring selection for one Q-ring"""

    def test_accepted(self):
        """Only seed and rescued rings contribute to the wedge angles"""
        accepted = {status for status in RingStatus if status.accepted}
        assert accepted == {RingStatus.SEED, RingStatus.RESCUED}


class TestPeakFit:
    """One Gaussian peak fitted to the azimuthal profile of a Q-ring"""

    @staticmethod
    def _peak(**changes) -> PeakFit:
        """Peak at 181 deg with height 0.2 +/- 0.002 and FWHM 17 +/- 0.5 deg, with optional field changes"""
        values = dict(height=0.2, height_error=0.002, center=181.0, center_error=0.1, fwhm=17.0, fwhm_error=0.5)
        values.update(changes)
        return PeakFit(**values)

    def test_from_fit_result(self):
        """Height, center and FWHM with their errors are read from the fit parameters of one Gaussian"""
        # Arrange: a negative Sigma describes the same peak, since the Gaussian depends only on Sigma squared
        result = {"f2.Height": (0.2, 0.002), "f2.PeakCentre": (361.3, 0.1), "f2.Sigma": (-7.5, 0.2)}
        # Act
        peak = PeakFit.from_fit_result(result, "f2")
        # Assert
        assert (peak.height, peak.height_error, peak.center, peak.center_error) == (0.2, 0.002, 361.3, 0.1)
        assert peak.fwhm == pytest.approx(7.5 * SIGMA_TO_FWHM)
        assert peak.fwhm_error == pytest.approx(0.2 * SIGMA_TO_FWHM)

    def test_significance_and_contrast(self):
        """Significance H/sigma_H and contrast H/A0 are NaN when the uncertainty or background is not positive"""
        peak = self._peak()
        assert peak.significance == pytest.approx(100.0)
        assert peak.contrast(0.25) == pytest.approx(0.8)
        for background in (0.0, -0.1, np.nan, np.inf):
            assert np.isnan(peak.contrast(background))
        for height_error in (0.0, -0.002, np.nan):
            assert np.isnan(self._peak(height_error=height_error).significance)

    @pytest.mark.parametrize(
        "changes", [{}, {"height_error": 0.0}, {"center_error": -0.1}, {"fwhm_error": np.nan}, {"center": np.inf}]
    )
    def test_is_valid(self, changes):
        """A peak is valid only if all values are finite and all uncertainties are positive"""
        assert self._peak(**changes).is_valid == (not changes)

    def test_weights_and_with_center(self):
        """Averaging weights are inverse uncertainties, and with_center returns a copy with another center"""
        peak = self._peak()
        moved = peak.with_center(1.3)
        assert peak.center_weight == pytest.approx(10.0)
        assert peak.fwhm_weight == pytest.approx(2.0)
        assert moved.center == 1.3 and peak.center == 181.0  # the original is unchanged
        assert moved.fwhm == peak.fwhm and moved.height == peak.height


class TestRingFit:
    """Fit of the azimuthal profile of one Q-ring and its selection outcome"""

    @staticmethod
    def _fit_result(centers=(1.0, 181.0)) -> dict:
        """Fit result as returned by _fitSpectrum: background 0.25, two peaks of height 0.2 and Sigma 7.5 deg"""
        result = {"fit_function": "name=FlatBackground", "chisq": (1.5, 0.0), "f0.A0": (0.25, 0.01)}
        for label, center in zip(("f1", "f2"), centers):
            result.update({f"{label}.Height": (0.2, 0.002), f"{label}.PeakCentre": (center, 0.1)})
            result[f"{label}.Sigma"] = (7.5, 0.2)
        return result

    def test_from_fit_result(self):
        """A fitted ring starts as rejected, holds a copy of the fit result and exposes its parameters"""
        result = self._fit_result()
        ring = RingFit.from_fit_result(3, 0.05, result)
        assert ring.is_fitted and ring.status is RingStatus.REJECTED and ring.reasons == []
        assert (ring.index, ring.q, ring.background, ring.chisq) == (3, 0.05, 0.25, 1.5)
        assert [peak.center for peak in ring.peaks] == [1.0, 181.0]
        ring.parameters["chisq"] = (9.0, 0.0)
        assert result["chisq"] == (1.5, 0.0)  # the ring holds a copy of the fit result

    def test_failed(self):
        """A ring that could not be fitted has the fit-failed status, its reason and no parameters"""
        ring = RingFit.failed(4, 0.06, "empty peak search window")
        assert not ring.is_fitted and ring.status is RingStatus.FIT_FAILED
        assert ring.reasons == ["empty peak search window"] and ring.parameters == {}

    def test_evaluated_peaks(self):
        """Both peaks are evaluated, except the mirrored peak in symmetric mode"""
        assert len(RingFit.from_fit_result(0, 0.05, self._fit_result()).evaluated_peaks) == 2
        symmetric = RingFit.from_fit_result(0, 0.05, self._fit_result(), symmetric=True)
        assert symmetric.evaluated_peaks == symmetric.peaks[:1]  # the mirrored peak is not evaluated

    @pytest.mark.parametrize(
        "centers, reference, expected",
        [
            ((1.0, 181.0), (0.0, 180.0), (1.0, 181.0)),  # same order as the reference
            ((181.0, 1.0), (0.0, 180.0), (1.0, 181.0)),  # swapped order
            ((2.0, 181.0), (180.0, 359.0), (181.0, 362.0)),  # swapped, and unwrapped across 360 deg
        ],
    )
    def test_aligned_to(self, centers, reference, expected):
        """Peaks are matched to the reference centers regardless of fit order and unwrapped near the reference"""
        ring = RingFit.from_fit_result(0, 0.05, self._fit_result(centers))
        aligned = ring.aligned_to(reference)
        assert [peak.center for peak in aligned] == pytest.approx(expected)
        assert [peak.center for peak in ring.peaks] == list(centers)  # the ring is unchanged

    def test_to_legacy_dict(self):
        """The legacy dictionary holds the fit parameters, or the joined reasons of a failed fit"""
        result = self._fit_result()
        assert RingFit.from_fit_result(0, 0.05, result).to_legacy_dict() == {**result, "error": None}
        failed = RingFit.failed(1, 0.06, "a")
        failed.reasons.append("b")
        assert failed.to_legacy_dict() == {"error": "a; b"}
