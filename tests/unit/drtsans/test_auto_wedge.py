# Unit tests for drtsans.auto_wedge
import numpy as np
import pytest

from drtsans.auto_wedge import (
    SIGMA_TO_FWHM,
    AcceptanceCriteria,
    ConsensusGroup,
    PeakFit,
    RingFit,
    RingStatus,
    _average_peaks,
    _circular_difference,
    _find_consensus,
)


@pytest.mark.parametrize(
    "angle, reference, expected",
    [(361.0, 1.0, 0.0), (1.0, 359.0, 2.0), (359.0, 1.0, -2.0), (181.0, 1.0, -180.0), (90.0, 450.0, 0.0)],
)
def test_circular_difference(angle, reference, expected):
    """Signed angle differences are wrapped into [-180, 180) deg"""
    assert _circular_difference(angle, reference) == pytest.approx(expected)


def _fit_result(centers=(1.0, 181.0), background=0.25, chisq=1.5, f1=None, f2=None) -> dict:
    """Fit result as returned by _fitSpectrum, with two peaks of height 0.2 +/- 0.002 and Sigma 7.5 +/- 0.2 deg.

    ``f1`` and ``f2`` optionally replace parameters of one peak, as in ``f1={"Height": (0.03, 0.005)}``.
    """
    result = {"fit_function": "name=FlatBackground", "chisq": (chisq, 0.0), "f0.A0": (background, 0.01)}
    for label, center, changes in zip(("f1", "f2"), centers, (f1, f2)):
        peak = {"Height": (0.2, 0.002), "PeakCentre": (center, 0.1), "Sigma": (7.5, 0.2), **(changes or {})}
        result.update({f"{label}.{name}": value for name, value in peak.items()})
    return result


def _ring(index, centers, center_error=0.1, fwhm=18.0, symmetric=False) -> RingFit:
    """Fitted Q-ring with the given peak centers, center uncertainty and FWHM, all in degrees"""
    peaks = [{"PeakCentre": (center, center_error), "Sigma": (fwhm / SIGMA_TO_FWHM, 0.2)} for center in centers]
    return RingFit.from_fit_result(
        index, 0.01 * (index + 1), _fit_result(centers, f1=peaks[0], f2=peaks[1]), symmetric
    )


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

    def test_from_fit_result(self):
        """A fitted ring starts as rejected, holds a copy of the fit result and exposes its parameters"""
        result = _fit_result()
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
        assert len(RingFit.from_fit_result(0, 0.05, _fit_result()).evaluated_peaks) == 2
        symmetric = RingFit.from_fit_result(0, 0.05, _fit_result(), symmetric=True)
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
        ring = RingFit.from_fit_result(0, 0.05, _fit_result(centers))
        aligned = ring.aligned_to(reference)
        assert [peak.center for peak in aligned] == pytest.approx(expected)
        assert [peak.center for peak in ring.peaks] == list(centers)  # the ring is unchanged

    def test_weight(self):
        """The weight is the sum of 1 / center_error over the evaluated peaks"""
        result = _fit_result(f1={"PeakCentre": (1.0, 0.1)}, f2={"PeakCentre": (181.0, 0.2)})
        assert RingFit.from_fit_result(0, 0.05, result).weight == pytest.approx(15.0)
        assert RingFit.from_fit_result(0, 0.05, result, symmetric=True).weight == pytest.approx(10.0)

    def test_to_legacy_dict(self):
        """The legacy dictionary holds the fit parameters, or the joined reasons of a failed fit"""
        result = _fit_result()
        assert RingFit.from_fit_result(0, 0.05, result).to_legacy_dict() == {**result, "error": None}
        failed = RingFit.failed(1, 0.06, "a")
        failed.reasons.append("b")
        assert failed.to_legacy_dict() == {"error": "a; b"}


class TestAcceptanceCriteria:
    """Thresholds deciding whether a Q-ring is a seed candidate or can be rescued"""

    @staticmethod
    def _criteria(**changes) -> AcceptanceCriteria:
        """Criteria with N = 10, C = 0.2, tolerance 0.5 FWHM, 1 deg azimuthal bins and maxchisq 1000"""
        settings = dict(signal_to_noise_min=10.0, contrast_min=0.2, center_tolerance=0.5, azimuthal_delta=1.0)
        settings.update(maxchisq=1000.0, **changes)
        return AcceptanceCriteria.from_settings(**settings)

    def test_from_settings(self):
        """The FWHM range spans 3 azimuthal bins to 180 deg, and the rescue threshold never exceeds N"""
        criteria = self._criteria(azimuthal_delta=2.0)
        assert (criteria.fwhm_min, criteria.fwhm_max, criteria.rescue_significance_min) == (6.0, 180.0, 5.0)
        assert criteria.consensus_weight_fraction == 0.67
        assert self._criteria(signal_to_noise_min=1.2).rescue_significance_min == 1.2

    @pytest.mark.parametrize(
        "changes",
        [
            {"contrast_min": 0.0},
            {"signal_to_noise_min": np.nan},
            {"center_tolerance": 0.0},
            {"center_tolerance": 1.5},
            {"azimuthal_delta": 60.0},
        ],
    )
    def test_invalid_settings(self, changes):
        """Non-positive thresholds, a tolerance outside (0, 1] and too wide azimuthal bins are refused"""
        with pytest.raises(ValueError):
            self._criteria(**changes)
        with pytest.raises(ValueError, match="Consensus weight fraction"):
            AcceptanceCriteria(10.0, 0.2, 0.5, 3.0, 1000.0, consensus_weight_fraction=0.5)

    def test_seed_failures_none(self):
        """A strong ring with H/sigma_H = 100, H/A0 = 0.8 and FWHM = 17.7 deg is a seed candidate"""
        assert self._criteria().seed_failures(RingFit.from_fit_result(0, 0.05, _fit_result())) == []

    @pytest.mark.parametrize(
        "changes, expected",
        [
            ({"chisq": 2000.0}, "chisq = 2e+03 > 1000"),
            ({"f1": {"Height": (0.2, 0.0)}}, "peak 1: non-finite parameter or non-positive uncertainty"),
            ({"f2": {"Height": (0.04, 0.002)}}, "peak 2: H/A0 = 0.16 < 0.2"),
            ({"background": -0.1}, "peak 1: background A0 = -0.1 is not positive"),
            ({"f1": {"Height": (0.2, 0.05)}}, "peak 1: H/sigma_H = 4 < 10"),
            ({"f2": {"Sigma": (0.5, 0.05)}}, "peak 2: FWHM = 1.18 deg outside [3, 180] deg"),
            ({"f1": {"Sigma": (80.0, 1.0)}}, "peak 1: FWHM = 188 deg outside [3, 180] deg"),
        ],
    )
    def test_seed_failures(self, changes, expected):
        """Each failed seed criterion is reported with the measured value and the threshold"""
        assert expected in self._criteria().seed_failures(RingFit.from_fit_result(0, 0.05, _fit_result(**changes)))

    def test_rescue_failures(self):
        """A weak ring of low contrast can be rescued if H/sigma_H >= 5, but not below"""
        criteria = self._criteria()
        weak = RingFit.from_fit_result(0, 0.05, _fit_result(f1={"Height": (0.035, 0.005)}))
        assert len(criteria.seed_failures(weak)) == 2  # contrast 0.14 and significance 7
        assert criteria.rescue_failures(weak) == []
        weaker = RingFit.from_fit_result(0, 0.05, _fit_result(f1={"Height": (0.015, 0.005)}))
        assert criteria.rescue_failures(weaker) == ["peak 1: H/sigma_H = 3 < 5"]

    def test_symmetric_and_failed_rings(self):
        """Only the fitted peak of a symmetric ring is checked, and a failed ring reports the failed fit"""
        criteria = self._criteria()
        symmetric = RingFit.from_fit_result(0, 0.05, _fit_result(f2={"Height": (0.01, 0.01)}), symmetric=True)
        assert criteria.seed_failures(symmetric) == []
        failed = RingFit.failed(1, 0.06, "empty peak search window")
        assert criteria.seed_failures(failed) == criteria.rescue_failures(failed) == ["fit failed"]


class TestConsensusGroup:
    """Seed candidates sharing one orientation, and the rules to match a ring to it"""

    GROUP = ConsensusGroup(centers=(181.0, 1.0), fwhms=(18.0, 18.0), members=(), weight_fraction=1.0)

    @pytest.mark.parametrize("centers, expected", [((186.0, 6.0), True), ((193.0, 13.0), False), ((6.0, 186.0), True)])
    def test_matches_rescue(self, centers, expected):
        """A rescued ring must lie within 0.5 x the group FWHM (9 deg) of the group center, in any peak order"""
        assert self.GROUP.matches(_ring(0, centers), tolerance=0.5) == expected

    def test_matches_candidate(self):
        """A seed candidate is compared using the mean of its FWHM and the group FWHM"""
        ring = _ring(0, (192.0, 12.0), fwhm=30.0)  # 11 deg away; tolerance 0.5 x (30 + 18) / 2 = 12 deg
        assert self.GROUP.matches(ring, tolerance=0.5, candidate=True)
        assert not self.GROUP.matches(ring, tolerance=0.5)  # rescue tolerance 0.5 x 18 = 9 deg

    def test_matches_symmetric(self):
        """Only the fitted peak of a symmetric ring is compared with the group center"""
        assert self.GROUP.matches(_ring(0, (186.0, 50.0), symmetric=True), tolerance=0.5)
        assert not self.GROUP.matches(_ring(0, (186.0, 50.0)), tolerance=0.5)


def test_average_peaks():
    """Centers and FWHMs are averaged with 1 / error weights, after aligning the peak order and the 360 deg wrap"""
    rings = [_ring(0, (180.0, 0.0), center_error=0.1), _ring(1, (183.0, 3.0), center_error=0.2)]
    expected_centers = pytest.approx((181.0, 1.0))  # weights 10 and 5
    assert _average_peaks(rings, (180.0, 0.0)) == (expected_centers, pytest.approx((18.0, 18.0)))
    rings[1] = _ring(1, (363.0, 183.0), center_error=0.2)  # swapped peak order, second center wrapped
    assert _average_peaks(rings, (180.0, 0.0))[0] == expected_centers
    with pytest.raises(ValueError):
        _average_peaks([], (180.0, 0.0))


def test_find_consensus_single_group():
    """Rings within 0.5 FWHM of the group center form one group, listed strongest first"""
    rings = [_ring(2, (186.0, 6.0), 0.25), _ring(1, (176.0, 356.0), 0.2), _ring(0, (181.0, 1.0), 0.1)]
    (group,) = _find_consensus(rings, tolerance=0.5)
    assert [ring.index for ring in group.members] == [0, 1, 2]
    assert group.weight_fraction == pytest.approx(1.0)


def test_find_consensus_two_groups():
    """Rings with different orientations form separate groups, largest weight fraction first"""
    rings = [_ring(0, (91.0, 271.0)), _ring(1, (181.0, 1.0)), _ring(2, (181.5, 1.5))]
    groups = _find_consensus(rings, tolerance=0.5)
    assert [[ring.index for ring in group.members] for group in groups] == [[1, 2], [0]]
    assert [group.weight_fraction for group in groups] == pytest.approx([2 / 3, 1 / 3])


def test_find_consensus_no_chaining():
    """Rings whose centers drift by 7 deg per ring are not chained into a single group"""
    rings = [_ring(index, (center, center + 180.0)) for index, center in enumerate((180.0, 187.0, 194.0, 201.0))]
    groups = _find_consensus(rings, tolerance=0.5)
    assert [[ring.index for ring in group.members] for group in groups] == [[0, 1], [2, 3]]
