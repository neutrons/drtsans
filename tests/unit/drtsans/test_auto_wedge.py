# Unit tests for drtsans.auto_wedge
import numpy as np
import pytest

import drtsans.auto_wedge
from drtsans.auto_wedge import (
    SIGMA_TO_FWHM,
    AcceptanceCriteria,
    AutoWedgeError,
    ConsensusGroup,
    PeakFit,
    RingFit,
    RingSelection,
    RingStatus,
    SelectionFailure,
    _average_peaks,
    _circular_difference,
    _effective_background_width,
    _failure_message,
    _find_consensus,
    _format_ring_table,
    _select_rings,
    _summarize_reasons,
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


def _ring(index, centers, center_error=0.1, fwhm=18.0, symmetric=False, heights=((0.2, 0.002),) * 2) -> RingFit:
    """Fitted Q-ring with the given peak centers, center uncertainty and FWHM in degrees, and peak heights"""
    peaks = [
        {"PeakCentre": (center, center_error), "Sigma": (fwhm / SIGMA_TO_FWHM, 0.2), "Height": height}
        for center, height in zip(centers, heights)
    ]
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


CRITERIA = AcceptanceCriteria.from_settings(10.0, 0.2, 0.5, 1.0, 1000.0)  # N, C, tolerance, bin width, maxchisq
WEAK = ((0.035, 0.005),) * 2  # H/A0 = 0.14 and H/sigma_H = 7: fails the seed criteria, passes the rescue criteria


def _statuses(selection: RingSelection) -> list[str]:
    return [ring.status.name for ring in selection.ring_fits]


class TestRingSelection:
    """Outcome of the Q-ring selection"""

    def test_properties(self):
        """Accepted rings are the seed and rescued rings, and legacy dictionaries are keyed by ring index"""
        rings = [_ring(0, (181.0, 1.0)), _ring(1, (182.0, 2.0)), RingFit.failed(2, 0.03, "no data in ring")]
        rings[0].status, rings[1].status = RingStatus.SEED, RingStatus.RESCUED
        selection = RingSelection(rings, CRITERIA, centers=(181.5, 1.5), fwhms=(18.0, 18.0))
        assert selection.succeeded and selection.accepted == rings[:2]
        legacy = selection.legacy_dicts()
        assert sorted(legacy) == [0, 1, 2]
        assert legacy[0]["error"] is None and legacy[2] == {"error": "no data in ring"}


def test_select_rings_no_fitted_rings():
    """Without any fitted ring the selection fails, and the rings keep their fit failure reasons"""
    selection = _select_rings([RingFit.failed(0, 0.01, "no data in ring")], CRITERIA)
    assert selection.failure is SelectionFailure.NO_FITTED_RINGS and selection.centers is None
    assert _statuses(selection) == ["FIT_FAILED"] and selection.ring_fits[0].reasons == ["no data in ring"]


def test_select_rings_no_seed_candidate():
    """Without any seed candidate the selection fails, and each ring reports its seed failures"""
    selection = _select_rings([_ring(0, (181.0, 1.0), heights=WEAK)], CRITERIA)
    assert selection.failure is SelectionFailure.NO_SEED_CANDIDATE and selection.centers is None
    assert _statuses(selection) == ["REJECTED"]
    assert "peak 1: H/sigma_H = 7 < 10" in selection.ring_fits[0].reasons


def test_select_rings_no_consensus():
    """Two equally strong orientations fail the 67 % rule, and both candidates are inconsistent"""
    selection = _select_rings([_ring(0, (181.0, 1.0)), _ring(1, (91.0, 271.0))], CRITERIA)
    assert selection.failure is SelectionFailure.NO_CONSENSUS and selection.centers is None
    assert [group.weight_fraction for group in selection.groups] == pytest.approx([0.5, 0.5])
    assert _statuses(selection) == ["INCONSISTENT", "INCONSISTENT"]


def _mixed_rings() -> list[RingFit]:
    """Three seeds near 181 deg, an inconsistent ring, a rescued ring, two rejected rings and a failed fit"""
    return [
        _ring(0, (181.0, 1.0)),
        _ring(1, (182.0, 2.0)),
        _ring(2, (180.0, 0.0)),
        _ring(3, (91.0, 271.0), center_error=0.3),
        _ring(4, (185.0, 5.0), center_error=1.0, heights=WEAK),
        _ring(5, (200.0, 20.0), center_error=1.0, heights=WEAK),
        _ring(6, (181.0, 1.0), heights=((0.01, 0.005),) * 2),
        RingFit.failed(7, 0.08, "empty peak search window"),
    ]


def test_select_rings_success():
    """Seeds, inconsistent, rescued and rejected rings are identified, and accepted rings are averaged"""
    selection = _select_rings(_mixed_rings(), CRITERIA)
    assert selection.succeeded
    expected = ["SEED"] * 3 + ["INCONSISTENT", "RESCUED", "REJECTED", "REJECTED", "FIT_FAILED"]
    assert _statuses(selection) == expected
    reasons = [ring.reasons for ring in selection.ring_fits]
    assert reasons[3] == ["orientation 91.0/271.0 deg differs from the consensus 181.0/1.0 deg"]
    assert reasons[5][-1] == "centers outside the consensus tolerance (181.0/1.0 deg)"
    assert reasons[6][-1] == "rescue: peak 2: H/sigma_H = 2 < 5"
    # weights 1 / center_error: 10 for each seed and 1 for the rescued ring at 185 deg
    assert selection.centers == pytest.approx((5615.0 / 31.0, 5615.0 / 31.0 - 180.0))
    assert selection.fwhms == pytest.approx((18.0, 18.0))


def test_select_rings_repeatable():
    """Selecting the same rings twice gives the same statuses and reasons"""
    rings = _mixed_rings()
    first = _select_rings(rings, CRITERIA)
    statuses, reasons = _statuses(first), [list(ring.reasons) for ring in first.ring_fits]
    second = _select_rings(rings, CRITERIA)
    assert _statuses(second) == statuses and [ring.reasons for ring in second.ring_fits] == reasons


def test_select_rings_symmetric():
    """In symmetric mode only the fitted peak is evaluated, so a weak mirror peak does not reject a ring"""
    heights = ((0.2, 0.002), (0.01, 0.01))
    rings = [_ring(index, (181.0, 1.0), symmetric=True, heights=heights) for index in range(2)]
    selection = _select_rings(rings, CRITERIA)
    assert selection.succeeded and _statuses(selection) == ["SEED", "SEED"]


class TestAutoWedgeError:
    """Error raised when automatic wedge detection fails"""

    def test_attributes(self):
        """The error is a RuntimeError, keeps the selection, and is exported by the module"""
        selection = _select_rings([RingFit.failed(0, 0.01, "no data in ring")], CRITERIA)
        error = AutoWedgeError("Automatic wedge detection failed", selection)
        assert isinstance(error, RuntimeError) and str(error) == "Automatic wedge detection failed"
        assert error.selection is selection
        assert "AutoWedgeError" in drtsans.auto_wedge.__all__


@pytest.mark.parametrize(
    "reasons, expected",
    [
        (["chisq = 2e+03 > 1000"], "chisq"),
        (["peak 1: non-finite parameter or non-positive uncertainty"], "invalid fit"),
        (["peak 1: background A0 = -0.1 is not positive"], "background"),
        (["peak 1: H/A0 = 0.1 < 0.2", "peak 2: H/A0 = 0.1 < 0.2"], "contrast"),  # named once for both peaks
        (["peak 2: H/sigma_H = 4 < 10"], "significance"),
        (["peak 1: FWHM = 1.5 deg outside [3, 180] deg"], "FWHM"),
        (["centers outside the consensus tolerance (181.0/1.0 deg)"], "center tolerance"),
        (["orientation 91.0/271.0 deg differs from the consensus 181.0/1.0 deg"], "orientation"),
        (["peak 1: H/sigma_H = 2 < 10", "rescue: peak 1: H/sigma_H = 2 < 5"], "significance; rescue: significance"),
    ],
)
def test_summarize_reasons(reasons, expected):
    """Each reason of a fitted ring is reduced to the name of its criterion"""
    ring = _ring(0, (181.0, 1.0))
    ring.reasons = reasons
    assert _summarize_reasons(ring) == expected


def test_summarize_reasons_fit_failed():
    """The reason of a failed fit is reported in full"""
    assert _summarize_reasons(RingFit.failed(0, 0.01, "Fit did not converge")) == "Fit did not converge"


def test_format_ring_table():
    """The table has a header and one row per ring in Q order, ending with the status and failed criteria"""
    rows = _format_ring_table(_select_rings(_mixed_rings(), CRITERIA)).splitlines()
    assert rows[0].split() == ["Q", "(1/A)", "H/A0", "H/sigma_H", "FWHM", "(deg)", "centers", "(deg)", "result"]
    results = ["seed"] * 3 + [
        "inconsistent: orientation",
        "rescued: contrast, significance",
        "rejected: contrast, significance, center tolerance",
        "rejected: contrast, significance; rescue: significance",
        "fit failed: empty peak search window",
    ]
    assert len(rows) == 9 and all(row.endswith(result) for row, result in zip(rows[1:], results))
    assert rows[1].split()[:8] == ["0.01", "0.8", "/", "0.8", "100", "/", "100", "18.0"]
    assert rows[-1].split()[:5] == ["0.08", "-", "-", "-", "-"]
    symmetric = _ring(0, (181.0, 1.0), symmetric=True)
    assert _format_ring_table(_select_rings([symmetric], CRITERIA)).splitlines()[1].split()[:5] == [
        "0.01",
        "0.8",
        "100",
        "18.0",
        "181.0",
    ]  # one value per column for the single fitted peak


@pytest.mark.parametrize(
    "rings, expected_lines",
    [
        (
            [_ring(0, (181.0, 1.0), heights=WEAK)],
            ["contrast H/A0 >= 0.2 (autoWedgeContrastMin)", "Highest contrast found: 0.14 at Q = 0.01 1/A"],
        ),
        (
            [_ring(0, (181.0, 1.0)), _ring(1, (91.0, 271.0))],
            ["Consensus tolerance: 0.5 x FWHM of each peak", "50%", "No group holds 67% of the seed weight."],
        ),
        ([RingFit.failed(0, 1.025, "no data in ring")], ["Adjust autoWedgeQmin / autoWedgeQmax"]),
    ],
)
def test_failure_message(rings, expected_lines):
    """Each kind of failure is explained with its key information and the per-ring table"""
    selection = _select_rings(rings, CRITERIA)
    message = _failure_message(selection)
    assert message.startswith(f"Automatic wedge detection failed: {selection.failure.value}.")
    assert all(line in message for line in expected_lines)
    assert _format_ring_table(selection) in message
    assert "Diagnostic plot" not in message


def test_failure_message_output_dir(tmp_path):
    """The paths of the diagnostic files are given when the output directory is known"""
    selection = _select_rings([RingFit.failed(0, 1.025, "no data in ring")], CRITERIA)
    message = _failure_message(selection, output_dir=str(tmp_path))
    assert f"{tmp_path / 'auto_wedge_rings.png'}, fit results: {tmp_path / 'auto_wedge_fit.h5'}" in message


@pytest.mark.parametrize(
    "centers, fwhms, expected",
    [
        ((180.0, 0.0), (18.0, 18.0), 1.5),  # narrow peaks: the requested width is kept
        ((180.0, 0.0), (59.9, 59.9), 1.5),  # just below the limit 1.5 x (FWHM1 + FWHM2) < 180 deg
        ((180.0, 0.0), (60.0, 60.0), 180.0 / 120.0 - 0.25),  # at the limit the width is reduced
        ((182.0, 3.0), (95.0, 94.7), 179.0 / 189.7 - 0.25),  # separation taken from the fitted centers
        ((358.0, 179.0), (95.0, 94.7), 179.0 / 189.7 - 0.25),  # order of the peaks and the 0/360 wrap
    ],
)
def test_effective_background_width(centers, fwhms, expected):
    """The background width is reduced only when the background wedges would be empty"""
    assert _effective_background_width(centers, fwhms, background_width=1.5, peak_width=0.25) == pytest.approx(
        expected
    )


def test_effective_background_width_no_room():
    """Peaks 180 deg wide leave no room for a background wedge outside the peak wedges"""
    assert _effective_background_width((180.0, 0.0), (180.0, 180.0), background_width=1.5, peak_width=0.25) <= 0.25
