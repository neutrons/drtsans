import numpy as np
import os
import re
import textwrap
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import List, Tuple
from drtsans.dataobjects import I1DAnnular
from drtsans.determine_bins import determine_1d_linear_bins
from drtsans.iq import BinningMethod, BinningParams, bin_annular_into_q1d

# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import DeleteWorkspace, logger, Gaussian, FlatBackground, mtd

# https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html
from mantid.simpleapi import Fit
from mantid.api import IFunction
import h5py
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

__all__ = ["AutoWedgeError", "getWedgeSelection"]

# factor to convert from a sigma to full width half max
# https://docs.mantidproject.org/nightly/fitting/fitfunctions/Gaussian.html
SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))

# Constants of the Q-ring selection
RESCUE_SIGNIFICANCE_MIN = 5.0  # H/sigma_H for a ring to confirm an orientation already found by the seed rings
FWHM_MIN_BINS = 3  # a peak narrower than three azimuthal bins is not resolved
FWHM_MAX = 180.0  # deg; the two peaks are 180 deg apart, so a wider peak is not defined
CONSENSUS_WEIGHT_FRACTION = 0.67  # fraction of the seed weight that the consensus group must hold


def _circular_difference(angle: float, reference: float) -> float:
    """Signed difference ``angle - reference`` between azimuthal angles, wrapped into [-180, 180) degrees.

    For example, 361 and 1 degrees differ by 0 degrees, and 1 and 359 degrees differ by 2 degrees.

    Parameters
    ----------
    angle: float
        Azimuthal angle, in degrees; any multiple of 360 deg is allowed.
    reference: float
        Azimuthal angle to compare with, in degrees; any multiple of 360 deg is allowed.

    Returns
    -------
    float
        Signed difference, in degrees, in [-180, 180).
    """
    return (angle - reference + 180.0) % 360.0 - 180.0


class RingStatus(Enum):
    """Outcome of the automatic wedge selection for one Q-ring."""

    FIT_FAILED = "fit failed"  # no usable fit of the azimuthal profile
    REJECTED = "rejected"  # fitted, but fails the seed criteria and is not rescued
    SEED = "seed"  # passes the seed criteria and belongs to the consensus orientation
    INCONSISTENT = "inconsistent"  # passes the seed criteria, but its orientation differs from the consensus
    RESCUED = "rescued"  # fails the seed criteria, but matches the consensus and passes the rescue criteria

    @property
    def accepted(self) -> bool:
        """True if rings with this status contribute to the wedge angles."""
        return self in (RingStatus.SEED, RingStatus.RESCUED)


@dataclass(frozen=True)
class PeakFit:
    """One Gaussian peak fitted to the azimuthal intensity profile I(phi) of a Q-ring.

    Attributes
    ----------
    height, height_error: float
        Peak height H above the flat background and its uncertainty, in the intensity units of the input
        (1/cm on absolute scale).
    center, center_error: float
        Peak center and its uncertainty, in degrees.
    fwhm, fwhm_error: float
        Full width at half maximum and its uncertainty, in degrees.
    """

    height: float
    height_error: float
    center: float
    center_error: float
    fwhm: float
    fwhm_error: float

    @classmethod
    def from_fit_result(cls, result: dict, label: str) -> "PeakFit":
        """Create from fit parameters stored as ``{"f1.Height": (value, error), ...}``.

        Parameters
        ----------
        result: dict
            Fit parameters as ``(value, error)`` tuples, keyed by Mantid parameter name.
        label: str
            Label of the Gaussian in the composite fit function, ``"f1"`` or ``"f2"``.
        """
        height, height_error = result[f"{label}.Height"]
        center, center_error = result[f"{label}.PeakCentre"]
        sigma, sigma_error = result[f"{label}.Sigma"]
        # the fit does not constrain the sign of Sigma
        fwhm, fwhm_error = abs(float(sigma)) * SIGMA_TO_FWHM, abs(float(sigma_error)) * SIGMA_TO_FWHM
        return cls(float(height), float(height_error), float(center), float(center_error), fwhm, fwhm_error)

    @property
    def is_valid(self) -> bool:
        """True if all values and uncertainties are finite and all uncertainties are positive."""
        errors = (self.height_error, self.center_error, self.fwhm_error)
        values = (self.height, self.center, self.fwhm) + errors
        return bool(np.all(np.isfinite(values)) and min(errors) > 0.0)

    @property
    def significance(self) -> float:
        """Peak height relative to its uncertainty, H / sigma_H; NaN if the uncertainty is not finite and positive."""
        if not (np.isfinite(self.height_error) and self.height_error > 0.0):
            return np.nan
        return self.height / self.height_error

    def contrast(self, background: float) -> float:
        """Peak height relative to the flat background A0 of the fit, H / A0; NaN unless A0 is finite and positive.

        Parameters
        ----------
        background: float
            Flat background A0 of the fit, in the intensity units of the input.

        Returns
        -------
        float
            Contrast H / A0, dimensionless, or NaN if A0 is not finite and positive.
        """
        if not (np.isfinite(background) and background > 0.0):
            return np.nan
        return self.height / background

    @property
    def center_weight(self) -> float:
        """Weight of the center when averaging over Q-rings, 1 / center_error."""
        return 1.0 / self.center_error

    @property
    def fwhm_weight(self) -> float:
        """Weight of the FWHM when averaging over Q-rings, 1 / fwhm_error."""
        return 1.0 / self.fwhm_error

    def with_center(self, center: float) -> "PeakFit":
        """Copy of this peak with another center; used to unwrap centers near a reference angle.

        Parameters
        ----------
        center: float
            New peak center, in degrees.

        Returns
        -------
        PeakFit
            Copy with the new center and an unchanged center error.
        """
        return replace(self, center=center)


@dataclass
class RingFit:
    """Fit of the azimuthal profile I(phi) of one Q-ring, and the outcome of the ring selection for it.

    Only ``status`` and ``reasons`` change after creation; the ring selection sets them.

    Attributes
    ----------
    index: int
        Index of the Q-ring.
    q: float
        Center of the Q-ring, in 1/Angstrom.
    parameters: dict
        Fit parameters as ``(value, error)`` tuples keyed by Mantid parameter name, plus ``"chisq"`` and the
        ``"fit_function"`` string. Empty when the fit failed.
    symmetric: bool
        True if the second peak is the mirror image of the first (``f2`` = ``f1`` shifted by 180 deg)
        rather than an independent fit.
    status: RingStatus
        Outcome of the selection. A fitted ring is ``REJECTED`` until the selection accepts it.
    reasons: list[str]
        Why the ring failed to fit or was not accepted.
    """

    index: int
    q: float
    parameters: dict
    symmetric: bool = False
    status: RingStatus = RingStatus.REJECTED
    reasons: list[str] = field(default_factory=list)

    @classmethod
    def from_fit_result(cls, index: int, q: float, result: dict, symmetric: bool = False) -> "RingFit":
        """Create from the dictionary returned by ``_fitSpectrum``.

        Parameters
        ----------
        index: int
            Index of the Q-ring.
        q: float
            Center of the Q-ring, in 1/Angstrom.
        result: dict
            Fit result of ``_fitSpectrum``: ``(value, error)`` tuples keyed by Mantid parameter name, plus
            ``"chisq"`` and ``"fit_function"``. It is copied.
        symmetric: bool
            True if ``f2`` is the mirror image of ``f1`` (``autoSymmetricWedges``).

        Returns
        -------
        RingFit
            Fitted ring with status ``REJECTED`` until the selection evaluates it.
        """
        return cls(index, q, dict(result), symmetric)

    @classmethod
    def failed(cls, index: int, q: float, reason: str) -> "RingFit":
        """Create for a Q-ring whose azimuthal profile could not be fitted.

        Parameters
        ----------
        index: int
            Index of the Q-ring.
        q: float
            Center of the Q-ring, in 1/Angstrom.
        reason: str
            Why the fit failed, such as the message of the exception raised by ``_fitSpectrum``.

        Returns
        -------
        RingFit
            Ring with no fit parameters and status ``FIT_FAILED``.
        """
        return cls(index, q, {}, status=RingStatus.FIT_FAILED, reasons=[reason])

    @property
    def is_fitted(self) -> bool:
        """True if fit parameters are available."""
        return self.status is not RingStatus.FIT_FAILED

    @property
    def background(self) -> float:
        """Flat background A0, in the intensity units of the input."""
        return float(self.parameters["f0.A0"][0])

    @property
    def chisq(self) -> float:
        """Cost function value of the fit."""
        return float(self.parameters["chisq"][0])

    @property
    def peaks(self) -> tuple[PeakFit, PeakFit]:
        """The two peaks, in the order of the fit function."""
        return PeakFit.from_fit_result(self.parameters, "f1"), PeakFit.from_fit_result(self.parameters, "f2")

    @property
    def evaluated_peaks(self) -> tuple[PeakFit, ...]:
        """Peaks subject to the acceptance criteria; a mirrored peak carries no independent information."""
        return self.peaks[:1] if self.symmetric else self.peaks

    @property
    def weight(self) -> float:
        """Weight of the ring in the consensus, the sum of 1 / center_error over the evaluated peaks."""
        return sum(peak.center_weight for peak in self.evaluated_peaks)

    def aligned_to(self, reference: tuple[float, float]) -> tuple[PeakFit, PeakFit]:
        """The two peaks matched to two reference centers, in degrees.

        The fit may return the peaks in either order, so the order with the smaller total angular distance
        to the reference is chosen. Each center is then expressed within 180 deg of its reference center,
        so that centers near 0 and 360 deg can be averaged.

        Parameters
        ----------
        reference: tuple[float, float]
            Reference centers of the two peaks, in degrees.

        Returns
        -------
        tuple[PeakFit, PeakFit]
            Peaks matched to the first and the second reference center, with centers within 180 deg of them.
        """
        first, second = self.peaks

        def distance(peaks: tuple[PeakFit, PeakFit]) -> float:
            return sum(abs(_circular_difference(p.center, r)) for p, r in zip(peaks, reference))

        ordered = min((first, second), (second, first), key=distance)
        return tuple(p.with_center(r + _circular_difference(p.center, r)) for p, r in zip(ordered, reference))

    def to_legacy_dict(self) -> dict:
        """Fit result in the dictionary format of the diagnostic plots and HDF5 file.

        ``"error"`` is None for a fitted ring, otherwise it holds the reasons of the fit failure.
        """
        error = None if self.is_fitted else "; ".join(self.reasons)
        return {**self.parameters, "error": error}


@dataclass(frozen=True)
class AcceptanceCriteria:
    """Thresholds deciding whether a Q-ring is a seed ring, or can be rescued once the orientation is known.

    A fitted ring whose peaks all pass the seed criteria (contrast, significance, FWHM range and chisq) is a seed
    candidate. Only seed candidates take part in finding the consensus orientation: candidates in the consensus
    group become seed rings, and candidates in any other group are inconsistent and not used. Rescuing applies
    afterwards to the rings that failed the seed criteria:
    such a ring is still accepted if its two peak centers match
    the consensus and its two peaks pass the weaker rescue criteria (significance, FWHM range and chisq,
    no contrast).
    Seed and rescued rings together determine the wedge angles.

    Attributes
    ----------
    significance_min: float
        N, minimum significance H/sigma_H of each peak of a seed ring.
    contrast_min: float
        C, minimum contrast H/A0 of each peak of a seed ring.
    center_tolerance: float
        Maximum deviation of each peak center from the consensus, as a fraction of that peak's FWHM.
    fwhm_min, fwhm_max: float
        Allowed range of the FWHM, in degrees.
    maxchisq: float
        Maximum cost function value of the fit.
    rescue_significance_min: float
        Minimum significance H/sigma_H of each peak of a rescued ring.
    consensus_weight_fraction: float
        Fraction of the seed weight that the consensus group must hold.
    """

    significance_min: float
    contrast_min: float
    center_tolerance: float
    fwhm_min: float
    maxchisq: float
    fwhm_max: float = FWHM_MAX
    rescue_significance_min: float = RESCUE_SIGNIFICANCE_MIN
    consensus_weight_fraction: float = CONSENSUS_WEIGHT_FRACTION

    @classmethod
    def from_settings(
        cls,
        signal_to_noise_min: float,
        contrast_min: float,
        center_tolerance: float,
        azimuthal_delta: float,
        maxchisq: float,
    ) -> "AcceptanceCriteria":
        """Create from the user settings and the azimuthal bin width.

        The rescue threshold never exceeds the seed threshold N, so that a seed-quality ring is never refused rescue.

        Parameters
        ----------
        signal_to_noise_min: float
            N, minimum significance H/sigma_H of each peak of a seed ring (``autoWedgeSignalToNoiseMin``).
        contrast_min: float
            C, minimum contrast H/A0 of each peak of a seed ring (``autoWedgeContrastMin``).
        center_tolerance: float
            Maximum deviation of each peak center from the consensus, as a fraction of the FWHM, in (0, 1]
            (``autoWedgeCenterTolerance``).
        azimuthal_delta: float
            Azimuthal bin width of the ring profiles, in degrees (``autoWedgeAzimuthalDelta``). The minimum FWHM
            is ``FWHM_MIN_BINS`` bins.
        maxchisq: float
            Maximum cost function value of the fit.

        Returns
        -------
        AcceptanceCriteria
            Criteria with the constants of the module for the other thresholds.

        Raises
        ------
        ValueError
            If a threshold is out of range (see ``__post_init__``).
        """
        rescue_significance_min = min(RESCUE_SIGNIFICANCE_MIN, signal_to_noise_min)
        fwhm_min = FWHM_MIN_BINS * azimuthal_delta
        return cls(
            signal_to_noise_min, contrast_min, center_tolerance, fwhm_min, maxchisq, FWHM_MAX, rescue_significance_min
        )

    def __post_init__(self):
        for name in ("significance_min", "contrast_min", "fwhm_min", "maxchisq", "rescue_significance_min"):
            value = getattr(self, name)
            if not (np.isfinite(value) and value > 0.0):
                raise ValueError(f"Auto-wedge {name} must be a positive number, got {value}")
        if not 0.0 < self.center_tolerance <= 1.0:
            raise ValueError(
                f"Auto-wedge center tolerance must be a fraction of the FWHM in (0, 1], got {self.center_tolerance}"
            )
        if self.fwhm_min >= self.fwhm_max:
            raise ValueError(
                f"Auto-wedge minimum FWHM ({FWHM_MIN_BINS} azimuthal bins = {self.fwhm_min:g} deg) must be below "
                f"{self.fwhm_max:g} deg; decrease autoWedgeAzimuthalDelta"
            )
        if not 0.5 < self.consensus_weight_fraction <= 1.0:
            raise ValueError(f"Consensus weight fraction must be in (0.5, 1], got {self.consensus_weight_fraction}")

    def seed_failures(self, ring: RingFit) -> list[str]:
        """Reasons why the ring is not a seed ring; empty if it is one.

        Parameters
        ----------
        ring: RingFit
            Q-ring to evaluate.

        Returns
        -------
        list[str]
            Failed seed criteria with the measured values, such as "peak 2: H/A0 = 0.16 < 0.2".
        """
        return self._failures(ring, self.significance_min, check_contrast=True)

    def rescue_failures(self, ring: RingFit) -> list[str]:
        """Reasons why the ring cannot be rescued, apart from its centers; empty if it can.

        Contrast is not required: the match with the consensus centers, checked separately, already protects
        against artifacts of low contrast.

        Parameters
        ----------
        ring: RingFit
            Q-ring to evaluate.

        Returns
        -------
        list[str]
            Failed rescue criteria with the measured values, such as "peak 1: H/sigma_H = 4 < 5".
        """
        return self._failures(ring, self.rescue_significance_min, check_contrast=False)

    def _failures(self, ring: RingFit, significance_min: float, check_contrast: bool) -> list[str]:
        """Failed criteria of a ring: chisq, then for each evaluated peak validity, contrast, significance and FWHM.

        Parameters
        ----------
        ring: RingFit
            Q-ring to evaluate.
        significance_min: float
            Minimum significance H/sigma_H of each peak.
        check_contrast: bool
            True to require a positive background and the minimum contrast ``contrast_min``.

        Returns
        -------
        list[str]
            Failed criteria with the measured values; ``["fit failed"]`` for a ring that was not fitted.
        """
        if not ring.is_fitted:
            return ["fit failed"]
        failures = [] if ring.chisq <= self.maxchisq else [f"chisq = {ring.chisq:.3g} > {self.maxchisq:g}"]
        for number, peak in enumerate(ring.evaluated_peaks, start=1):
            if not peak.is_valid:
                failures.append(f"peak {number}: non-finite parameter or non-positive uncertainty")
                continue
            contrast = peak.contrast(ring.background)
            if check_contrast and np.isnan(contrast):
                failures.append(f"peak {number}: background A0 = {ring.background:.3g} is not positive")
            elif check_contrast and contrast < self.contrast_min:
                failures.append(f"peak {number}: H/A0 = {contrast:.3g} < {self.contrast_min:g}")
            if peak.significance < significance_min:
                failures.append(f"peak {number}: H/sigma_H = {peak.significance:.3g} < {significance_min:g}")
            if not self.fwhm_min <= peak.fwhm <= self.fwhm_max:
                failures.append(
                    f"peak {number}: FWHM = {peak.fwhm:.3g} deg outside [{self.fwhm_min:g}, {self.fwhm_max:g}] deg"
                )
        return failures


@dataclass(frozen=True)
class ConsensusGroup:
    """Seed candidates sharing one orientation, with the weighted average of their peaks.

    Attributes
    ----------
    centers, fwhms: tuple[float, float]
        Weighted average center and FWHM of each of the two peaks, in degrees.
    members: tuple[RingFit, ...]
        Seed candidates in the group, strongest first.
    weight_fraction: float
        Share of the group in the total weight of all seed candidates.
    """

    centers: tuple[float, float]
    fwhms: tuple[float, float]
    members: tuple[RingFit, ...]
    weight_fraction: float

    def matches(self, ring: RingFit, tolerance: float, candidate: bool = False) -> bool:
        """True if each evaluated peak center of the ring lies within ``tolerance`` x FWHM of the group center.

        ``tolerance`` is a fraction of the FWHM of each peak. The FWHM is the group FWHM when rescuing a ring, and
        the mean of the ring and group FWHMs when a seed candidate (``candidate=True``) joins the group.

        Parameters
        ----------
        ring: RingFit
            Fitted Q-ring to compare with the group.
        tolerance: float
            Maximum deviation of each peak center, as a fraction of the FWHM (``autoWedgeCenterTolerance``).
        candidate: bool
            True when a seed candidate joins the group, False when a ring is rescued.

        Returns
        -------
        bool
            True if every evaluated peak center of the ring is within the tolerance.
        """
        aligned = ring.aligned_to(self.centers)[: len(ring.evaluated_peaks)]
        for peak, center, fwhm in zip(aligned, self.centers, self.fwhms):
            width = 0.5 * (peak.fwhm + fwhm) if candidate else fwhm
            if abs(peak.center - center) > tolerance * width:
                return False
        return True


def _average_peaks(
    rings: list[RingFit], reference: tuple[float, float]
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Weighted average centers and FWHMs of the two peaks over Q-rings, in degrees.

    The peaks of each ring are matched to the reference centers and unwrapped near them (see
    ``RingFit.aligned_to``). Centers are weighted by 1 / center_error and FWHMs by 1 / fwhm_error.

    Parameters
    ----------
    rings: list[RingFit]
        Fitted Q-rings to average, such as the accepted rings or the members of a consensus group.
    reference: tuple[float, float]
        Approximate centers of the two peaks, in degrees, to which the peaks of every ring are matched.

    Returns
    -------
    tuple[tuple[float, float], tuple[float, float]]
        ``(centers, fwhms)``: the weighted average center and FWHM of each of the two peaks, in degrees.
        The centers lie within 180 deg of the reference centers.

    Raises
    ------
    ValueError
        If ``rings`` is empty.
    """
    if not rings:
        raise ValueError("Cannot average the peaks of zero Q-rings")
    aligned = [ring.aligned_to(reference) for ring in rings]
    centers, fwhms = [], []
    for k in range(2):
        peaks = [pair[k] for pair in aligned]
        centers.append(float(np.average([p.center for p in peaks], weights=[p.center_weight for p in peaks])))
        fwhms.append(float(np.average([p.fwhm for p in peaks], weights=[p.fwhm_weight for p in peaks])))
    return (centers[0], centers[1]), (fwhms[0], fwhms[1])


def _find_consensus(candidates: list[RingFit], tolerance: float) -> list[ConsensusGroup]:
    """Group seed candidates by orientation, largest weight first.

    Candidates are taken strongest first. Each joins the first group whose current weighted center it matches
    (see ``ConsensusGroup.matches``), otherwise it starts a new group. The group center is recomputed after each
    addition. Since the strongest rings join first, they fix the center, which prevents chaining rings whose
    centers drift gradually with Q into one group.

    Parameters
    ----------
    candidates: list[RingFit]
        Seed candidates, the fitted Q-rings passing the seed criteria. Must not be empty.
    tolerance: float
        Maximum deviation of each peak center from the group center, as a fraction of the FWHM
        (``autoWedgeCenterTolerance``).

    Returns
    -------
    list[ConsensusGroup]
        Groups of candidates with the same orientation, largest weight fraction first.
    """
    groups: list[ConsensusGroup] = []
    for ring in sorted(candidates, key=lambda candidate: candidate.weight, reverse=True):
        index = next((i for i, g in enumerate(groups) if g.matches(ring, tolerance, candidate=True)), None)
        if index is None:
            centers = tuple(peak.center for peak in ring.peaks)
            groups.append(ConsensusGroup(centers, tuple(peak.fwhm for peak in ring.peaks), (ring,), 0.0))
        else:
            members = groups[index].members + (ring,)
            centers, fwhms = _average_peaks(list(members), groups[index].centers)
            groups[index] = ConsensusGroup(centers, fwhms, members, 0.0)
    total_weight = sum(ring.weight for ring in candidates)
    groups = [replace(g, weight_fraction=sum(ring.weight for ring in g.members) / total_weight) for g in groups]
    return sorted(groups, key=lambda group: group.weight_fraction, reverse=True)


class SelectionFailure(Enum):
    """Reason why no wedge angles could be determined from the Q-rings."""

    NO_FITTED_RINGS = "none of the Q-rings could be fitted"
    NO_SEED_CANDIDATE = "no Q-ring passes the seed criteria"
    NO_CONSENSUS = "seed candidates do not agree on a single orientation"


@dataclass
class RingSelection:
    """Outcome of the Q-ring selection: the status of every ring and, on success, the wedge peak parameters.

    Attributes
    ----------
    ring_fits: list[RingFit]
        All Q-rings, with their final status and reasons.
    criteria: AcceptanceCriteria
        Thresholds used for the selection.
    groups: list[ConsensusGroup]
        Groups of seed candidates by orientation, largest weight fraction first; the first is the consensus.
    centers, fwhms: tuple[float, float] or None
        Weighted average center and FWHM of the two peaks over the accepted rings, in degrees; None on failure.
    failure: SelectionFailure or None
        Reason of the failure; None on success.
    """

    ring_fits: list[RingFit]
    criteria: AcceptanceCriteria
    groups: list[ConsensusGroup] = field(default_factory=list)
    centers: tuple[float, float] | None = None
    fwhms: tuple[float, float] | None = None
    failure: SelectionFailure | None = None

    @property
    def succeeded(self) -> bool:
        """True if the wedge peak parameters were determined."""
        return self.failure is None

    @property
    def accepted(self) -> list[RingFit]:
        """Seed and rescued rings, which determine the wedge angles."""
        return [ring for ring in self.ring_fits if ring.status.accepted]

    def legacy_dicts(self) -> dict[int, dict]:
        """Fit results keyed by ring index, in the format of the diagnostic plots and HDF5 file."""
        return {ring.index: ring.to_legacy_dict() for ring in self.ring_fits}


def _select_rings(ring_fits: list[RingFit], criteria: AcceptanceCriteria) -> RingSelection:
    """Select the Q-rings that determine the wedge angles, setting the status and reasons of every ring.

    Seed candidates are the fitted rings passing the seed criteria. The consensus is the group of candidates
    holding at least ``criteria.consensus_weight_fraction`` of their weight; its members become seed rings. Rings
    failing the seed criteria are rescued if they match the consensus and pass the rescue criteria. A failed
    selection is reported through ``RingSelection.failure``, never raised.

    Parameters
    ----------
    ring_fits: list[RingFit]
        All Q-rings, fitted or not. The ``status`` and ``reasons`` of every ring are set in place.
    criteria: AcceptanceCriteria
        Thresholds of the seed and rescue criteria and of the consensus.

    Returns
    -------
    RingSelection
        Status and reasons of every Q-ring and, on success, the weighted average peak centers and FWHMs
        of the accepted rings.
    """
    selection = RingSelection(ring_fits, criteria)
    # Evaluate the seed criteria for every fitted ring, starting from a rejected status
    fitted = [ring for ring in ring_fits if ring.is_fitted]
    for ring in fitted:
        ring.status, ring.reasons = RingStatus.REJECTED, criteria.seed_failures(ring)
    if not fitted:
        selection.failure = SelectionFailure.NO_FITTED_RINGS
        return selection
    # Rings with no seed failures are the seed candidates
    candidates = [ring for ring in fitted if not ring.reasons]
    if not candidates:
        selection.failure = SelectionFailure.NO_SEED_CANDIDATE
        return selection

    # Group the candidates by orientation; the largest group must hold enough weight to be the consensus
    selection.groups = _find_consensus(candidates, criteria.center_tolerance)
    consensus = selection.groups[0]
    if consensus.weight_fraction < criteria.consensus_weight_fraction:
        for ring in candidates:
            ring.status = RingStatus.INCONSISTENT
            ring.reasons = [f"no orientation group holds {criteria.consensus_weight_fraction:.0%} of the seed weight"]
        selection.failure = SelectionFailure.NO_CONSENSUS
        return selection

    # Candidates in the consensus become seed rings, and those in other groups are inconsistent
    consensus_centers = "/".join(f"{center % 360.0:.1f}" for center in consensus.centers)
    for group in selection.groups:
        for ring in group.members:
            ring.status = RingStatus.SEED if group is consensus else RingStatus.INCONSISTENT
            if group is not consensus:
                centers = "/".join(f"{center % 360.0:.1f}" for center in group.centers)
                ring.reasons = [f"orientation {centers} deg differs from the consensus {consensus_centers} deg"]

    # Rescue rejected rings that match the consensus centers and pass the rescue criteria
    for ring in fitted:
        if ring.status is not RingStatus.REJECTED:
            continue
        rescue_failures = criteria.rescue_failures(ring)
        if not consensus.matches(ring, criteria.center_tolerance):
            ring.reasons.append(f"centers outside the consensus tolerance ({consensus_centers} deg)")
        elif rescue_failures:
            ring.reasons.extend(f"rescue: {reason}" for reason in rescue_failures)
        else:
            ring.status = RingStatus.RESCUED

    # Average the peaks of the seed and rescued rings, which determine the wedge angles
    selection.centers, selection.fwhms = _average_peaks(selection.accepted, consensus.centers)
    return selection


class AutoWedgeError(RuntimeError):
    """Automatic wedge detection could not determine the wedge angles.

    Raised when the Q-ring selection fails, or when the accepted peaks are too wide to leave room for background
    wedges; in the latter case ``selection.succeeded`` is True.

    Parameters
    ----------
    message: str
        Explanation of the failure, with the per-ring table and suggestions to the user.
    selection: RingSelection
        Outcome of the Q-ring selection.

    Attributes
    ----------
    selection: RingSelection
        Outcome of the Q-ring selection, with the status and reasons of every ring.
    """

    def __init__(self, message: str, selection: RingSelection):
        super().__init__(message)
        self.selection = selection


# Short names of the criteria, keyed by a text identifying the reasons produced by the ring selection
_CRITERION_NAMES = {
    "chisq": "chisq",
    "non-finite": "invalid fit",
    "background A0": "background",
    "H/A0": "contrast",
    "H/sigma_H": "significance",
    "FWHM": "FWHM",
    "centers outside": "center tolerance",
    "orientation": "orientation",
}


def _summarize_reasons(ring: RingFit) -> str:
    """Names of the criteria a ring failed, such as "contrast, significance; rescue: significance".

    The measured values are in the other table columns. A fit failure is reported in full.

    Parameters
    ----------
    ring: RingFit
        Q-ring after the selection, with its ``reasons`` set.

    Returns
    -------
    str
        Short names of the failed seed criteria, followed by those of the rescue criteria after "rescue:".
        Empty if the ring failed no criterion.
    """
    if not ring.is_fitted:
        return "; ".join(ring.reasons)
    seed, rescue = [], []
    for reason in ring.reasons:
        names = rescue if reason.startswith("rescue:") else seed
        name = next((n for key, n in _CRITERION_NAMES.items() if key in reason), reason)
        if name not in names:
            names.append(name)
    return ", ".join(seed) + ("; rescue: " + ", ".join(rescue) if rescue else "")


def _format_ring_table(selection: RingSelection) -> str:
    """Table of the fitted values and the selection outcome of every Q-ring, one row per ring in Q order.

    Parameters
    ----------
    selection: RingSelection
        Outcome of the Q-ring selection.

    Returns
    -------
    str
        Table with the Q-ring center (1/A), contrast H/A0, significance H/sigma_H, FWHM (deg) and center (deg)
        of each evaluated peak, and the ring status with the names of the failed criteria.
    """

    def per_peak(fmt: str, values: list[float]) -> str:
        return " / ".join(fmt.format(value) for value in values)

    table = [["Q (1/A)", "H/A0", "H/sigma_H", "FWHM (deg)", "centers (deg)", "result"]]
    for ring in sorted(selection.ring_fits, key=lambda ring_fit: ring_fit.q):
        columns = ["-"] * 4
        if ring.is_fitted:
            peaks = ring.evaluated_peaks
            columns = [
                per_peak("{:.3g}", [peak.contrast(ring.background) for peak in peaks]),
                per_peak("{:.3g}", [peak.significance for peak in peaks]),
                per_peak("{:.1f}", [peak.fwhm for peak in peaks]),
                per_peak("{:.1f}", [peak.center % 360.0 for peak in peaks]),
            ]
        result = ring.status.value + (": " + _summarize_reasons(ring) if ring.reasons else "")
        table.append([f"{ring.q:.4g}", *columns, result])
    # each column is as wide as its longest entry plus a two-space gap, so that wide values never run together
    widths = [max(len(row[column]) for row in table) + 2 for column in range(len(table[0]) - 1)]
    return "\n".join("   " + "".join(cell.ljust(width) for cell, width in zip(row, widths)) + row[-1] for row in table)


def _failure_message(selection: RingSelection, output_dir: str | None = None) -> str:
    """Explanation of a failed Q-ring selection, with the per-ring table and suggestions to the user.

    Parameters
    ----------
    selection: RingSelection
        Outcome of a failed Q-ring selection; ``selection.failure`` is not None.
    output_dir: str or None
        Directory of the diagnostic files, whose paths are added to the message. None leaves the paths out.

    Returns
    -------
    str
        Multi-line message, starting with "Automatic wedge detection failed:".
    """
    criteria = selection.criteria
    q_values = [ring.q for ring in selection.ring_fits]
    q_range = f"Q-ring centers {min(q_values):.4g}-{max(q_values):.4g} 1/A" if q_values else "no Q-rings"
    lines = [f"Automatic wedge detection failed: {selection.failure.value}."]
    if selection.failure is SelectionFailure.NO_SEED_CANDIDATE:
        lines += [
            f"  Seed criteria (each peak): contrast H/A0 >= {criteria.contrast_min:g} (autoWedgeContrastMin),",
            f"    significance H/sigma_H >= {criteria.significance_min:g} (autoWedgeSignalToNoiseMin),",
            f"    {criteria.fwhm_min:g} deg <= FWHM <= {criteria.fwhm_max:g} deg,",
            f"    chi2 <= {criteria.maxchisq:g} (each ring)",
            f"  {q_range}, {len(q_values)} rings (autoWedgeQmin, autoWedgeQmax, autoWedgeQdelta)",
            "",
            _format_ring_table(selection),
            "",
        ]
        contrasts = [
            (peak.contrast(ring.background), ring.q)
            for ring in selection.ring_fits
            if ring.is_fitted
            for peak in ring.evaluated_peaks
            if peak.is_valid and np.isfinite(peak.contrast(ring.background))
        ]
        if contrasts:
            contrast, q = max(contrasts)
            lines.append(
                f"  Highest contrast found: {contrast:.3g} at Q = {q:.4g} 1/A (required {criteria.contrast_min:g})."
            )
        lines += [
            "  Suggestions:",
            '   - Inspect I(Qx,Qy): the data may be isotropic in this Q range; use 1DQbinType "scalar"',
            "     or set WedgeMinAngles/WedgeMaxAngles manually.",
            "   - If the anisotropy lies at lower Q, decrease autoWedgeQmin.",
            "   - Lower autoWedgeContrastMin only if the anisotropy is known to be weak.",
        ]
    elif selection.failure is SelectionFailure.NO_CONSENSUS:
        lines += [
            f"  Consensus tolerance: {criteria.center_tolerance:g} x FWHM of each peak (autoWedgeCenterTolerance)",
            "",
            f"   {'group':<7}{'rings (Q, 1/A)':<30}{'centers (deg)':<16}weight",
        ]
        for label, group in zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", selection.groups):
            rings = ", ".join(f"{q:.4g}" for q in sorted(ring.q for ring in group.members))
            centers = " / ".join(f"{center % 360.0:.1f}" for center in group.centers)
            lines.append(f"   {label:<7}{rings:<30}{centers:<16}{group.weight_fraction:.0%}")
        lines += [
            "",
            f"  No group holds {criteria.consensus_weight_fraction:.0%} of the seed weight.",
            "  The sample may contain more than one oriented structure. Restrict the Q range to one",
            "  group (autoWedgeQmin / autoWedgeQmax), or set WedgeMinAngles/WedgeMaxAngles manually.",
            "",
            _format_ring_table(selection),
        ]
    else:
        lines += [
            f"  {q_range}, {len(q_values)} rings",
            "",
            _format_ring_table(selection),
            "",
            "  Adjust autoWedgeQmin / autoWedgeQmax to lie within the Q range of the data.",
        ]
    if output_dir is not None:
        plots, fits = os.path.join(output_dir, "auto_wedge_rings.png"), os.path.join(output_dir, "auto_wedge_fit.h5")
        lines.append(f"  Diagnostic plot of the ring fits: {plots}, fit results: {fits}")
    return "\n".join(lines)


def _effective_background_width(
    centers: tuple[float, float], fwhms: tuple[float, float], background_width: float, peak_width: float
) -> float:
    """Background width, as a fraction of the FWHM, leaving background wedges at least as wide as the peak wedges.

    Each background wedge spans from ``background_width`` x FWHM past one peak to ``background_width`` x FWHM
    before the other, so its width is S - ``background_width`` x (FWHM1 + FWHM2), with S the separation of the peak
    centers. The requested width is kept if the background wedges are at least as wide as the peak wedges, i.e.
    (``background_width`` + ``peak_width``) x (FWHM1 + FWHM2) <= S. Otherwise the width is reduced to
    S / (FWHM1 + FWHM2) - ``peak_width``, for which each background wedge is as wide as the peak wedges.

    Parameters
    ----------
    centers: tuple[float, float]
        Centers of the two peaks, in degrees. Their order and any multiple of 360 deg do not matter.
    fwhms: tuple[float, float]
        FWHM of the two peaks, in degrees.
    background_width: float
        Requested distance of each background wedge from the peak centers, as a fraction of the FWHM
        (``autoWedgeBackgroundWidth``).
    peak_width: float
        Half-width of each peak wedge, as a fraction of the FWHM (``autoWedgePeakWidth``).

    Returns
    -------
    float
        Background width to use, as a fraction of the FWHM. A reduced value not above ``peak_width`` means that the
        peaks leave no room for a background wedge.
    """
    difference = abs(_circular_difference(centers[1], centers[0]))
    separation = min(difference, 360.0 - difference)
    total_fwhm = fwhms[0] + fwhms[1]
    if (background_width + peak_width) * total_fwhm <= separation:
        return background_width
    return separation / total_fwhm - peak_width


def _create_fit_function(combo_function_str):
    """Create fit function from a Mantid combo function string

    An example:
    name=FlatBackground,A0=5894.8246694221725;name=Gaussian,Height=16626.818784787953,
    PeakCentre=100.43905821664455,Sigma=5.917594016360494;name=Gaussian,Height=18517.052761733517,
    PeakCentre=304.6715713990027,Sigma=36.943377164429435

    Parameters
    ----------
    combo_function_str: str
        function string

    Returns
    -------
    ~list
        list of Mantid IFunction in order of

    """
    mantid_function_map = {Gaussian().name: Gaussian, FlatBackground().name: FlatBackground}

    # split to multiple functions
    single_function_str_list = combo_function_str.split(";")
    single_functions = list()

    # convert to function
    for single_func in single_function_str_list:
        # split with ',' except inside parentheses, as in "constraints=(0<Height,0<Sigma)"
        params = re.split(r",(?![^()]*\))", single_func)
        # create function
        func_name = params[0].split("=")[1]
        func_i = mantid_function_map[func_name]()
        # start to set parameters
        for value_index in range(1, len(params)):
            # split in PARAM_NAME=PARAM_VALUE style
            param_list = params[value_index].split("=")
            param_name_i = param_list[0]
            param_value_i = param_list[1]
            if param_name_i == "constraints":
                for constraint in param_value_i.strip("()").split(","):
                    func_i.constrain(constraint)
            else:
                func_i.__setattr__(param_name_i, float(param_value_i))
        # END-FOR
        single_functions.append(func_i)

    return single_functions


def _set_function_param_value(functions, param_name, param_value):
    """set function parameter value from the result of a combo function

    Parameters
    ----------
    functions: ~list
        list of functions
    param_name: str
        function parameter name may have 2 forms: PARAM_NAME or fI.PARAM_NAME, where I is an integer
    param_value: float
        parameter value

    Returns
    -------
    None

    """
    # Check parameter name
    if "." in param_name:
        terms = param_name.split(".")
        function_index = int(terms[0].split("f")[-1])
        param_name = terms[1]
    else:
        function_index = 0

    # set
    functions[function_index].__setattr__(param_name, param_value)


def _calculate_function(fit_functions, vec_x):
    """Calculate a function (set)

    Parameters
    ----------
    fit_functions: ~list[IFunction]
        List of fit function
    vec_x: numpy.ndarray
        1D array for X values

    Returns
    -------
    numpy.ndarray
        1D array for calculated function values

    """
    # Create the combo function from a list of functions
    combo_function = fit_functions[0]
    for f_i in range(1, len(fit_functions)):
        combo_function += fit_functions[f_i]

    # Call the create
    vec_y = combo_function.__call__(vec_x)

    return vec_y


# Layout of the figure of the Q-ring fits
RING_PLOT_COLUMNS = 5
RING_PLOT_PANEL_SIZE = (4.25, 5.0)  # inches


def _ring_title(ring_fit: RingFit) -> str:
    """Q, status and failed criteria of a Q-ring, such as "Q = 0.01 1/A" and "rejected: contrast" on two lines.

    Parameters
    ----------
    ring_fit: RingFit
        Q-ring after the selection.

    Returns
    -------
    str
        Title of the ring's panel in ``auto_wedge_rings.png``. The second line is wrapped to 45 characters.
    """
    outcome = ring_fit.status.value + (f": {_summarize_reasons(ring_fit)}" if ring_fit.reasons else "")
    return f"Q = {ring_fit.q:.4g} 1/A\n" + textwrap.fill(outcome, width=45)


def _plot_fit_results(rings, selection: RingSelection, output_dir, auto_symmetric_wedges):
    """Plot the observed profile and the fit of every Q-ring in one figure, ``auto_wedge_rings.png``

    There is one panel per Q-ring in Q order, including rings whose fit failed (observed points only).
    The panel title gives the ring status and the criteria it failed.

    Parameters
    ----------
    rings:  ~list
        List of I(Qx, Qy) in rings
    selection: RingSelection
        Fit results and selection outcome of every Q-ring
    output_dir: str
        full path of the output directory
    auto_symmetric_wedges: bool
        if True, label the second peak as mirrored
    """

    def _set_function_params_and_calculate(fit_functions: list[IFunction], peak_fitted_dict: dict):
        for param_name in peak_fitted_dict:
            if param_name not in ["fit_function", "error", "used"]:
                # parameter value is recorded as tuple as value and error
                _set_function_param_value(fit_functions, param_name, peak_fitted_dict[param_name][0])
        return _calculate_function(fit_functions, ring.phi)

    peak_fit_dict = selection.legacy_dicts()
    ring_fits = {ring_fit.index: ring_fit for ring_fit in selection.ring_fits}
    rows = max(1, -(-len(rings) // RING_PLOT_COLUMNS))  # ceiling division
    width, height = RING_PLOT_PANEL_SIZE
    figure, axes = plt.subplots(rows, RING_PLOT_COLUMNS, figsize=(RING_PLOT_COLUMNS * width, rows * height))
    axes = np.atleast_1d(axes).ravel()
    for axis in axes[len(rings) :]:
        axis.set_visible(False)

    for index, (ring, axis) in enumerate(zip(rings, axes)):
        # markers, since a line is not drawn between bins separated by NaN (bins without data)
        axis.plot(ring.phi, ring.intensity, "o", markersize=3, label="observed", color="black", zorder=3)
        axis.set_title(_ring_title(ring_fits[index]), fontsize="large")
        axis.tick_params(labelsize="medium")
        axis.set_xlabel("azimuthal angle (deg)", fontsize="large")
        if index % RING_PLOT_COLUMNS == 0:
            axis.set_ylabel("intensity", fontsize="large")  # first panel of each row
        if peak_fit_dict[index]["error"] is not None:
            continue  # the fit failed; show the observed points only

        # construct the fit function; its parameters are set to the fitted values below
        fit_function_set = _create_fit_function(peak_fit_dict[index]["fit_function"])

        if auto_symmetric_wedges:
            # separate fitted peak ("f1" params) from mirrored peak (fake "f2" params) to plot them separately
            # the fit function is a single gaussian with only "f1" entries
            peak_fitted_dict = {}
            peak_mirrored_dict = {}
            for key, val in peak_fit_dict[index].items():
                if key.startswith("f2"):
                    new_key = key.replace("f2", "f1")
                    # shift mirrored peak into plot angle range
                    new_value = (val[0] - 360.0, val[1]) if "PeakCentre" in key and val[0] > 360.0 else val
                    peak_mirrored_dict[new_key] = new_value
                else:
                    peak_fitted_dict[key] = val
        else:
            peak_fitted_dict = peak_fit_dict[index]

        model_y = _set_function_params_and_calculate(fit_function_set, peak_fitted_dict)
        axis.plot(ring.phi, model_y, label="fitted", color="red")
        if auto_symmetric_wedges:
            mirrored_model_y = _set_function_params_and_calculate(fit_function_set, peak_mirrored_dict)
            axis.plot(ring.phi, mirrored_model_y, label="mirrored", color="red", linestyle="dashed")

    # one legend for the figure, taken from the panel with the most curves
    handles, labels = max((axis.get_legend_handles_labels() for axis in axes), key=lambda hl: len(hl[1]))
    figure.legend(handles, labels, loc="upper center", ncol=len(labels), fontsize="large", frameon=False)
    # leave 0.4 inch for the legend, and 3 font heights between rows for the separating line
    figure.tight_layout(rect=(0, 0, 1, 1 - 0.4 / (rows * height)), h_pad=3.0)
    # separate the rows with a horizontal line, in the middle of the gap between the x-axis labels of a row
    # and the titles of the next row
    to_figure = figure.transFigure.inverted()

    def _extent(row: int) -> list:
        visible = [
            axis for axis in axes[row * RING_PLOT_COLUMNS : (row + 1) * RING_PLOT_COLUMNS] if axis.get_visible()
        ]
        return [axis.get_tightbbox().transformed(to_figure) for axis in visible]

    for row in range(1, rows):
        y = 0.5 * (min(box.y0 for box in _extent(row - 1)) + max(box.y1 for box in _extent(row)))
        figure.add_artist(Line2D([0.01, 0.99], [y, y], transform=figure.transFigure, color="gray", linewidth=1))
    figure.savefig(os.path.join(output_dir, "auto_wedge_rings.png"))
    plt.close(figure)


def _export_to_h5(iq2d, rings, azimuthal_delta, selection: RingSelection, output_dir):
    """write annular binned I(Qx, Qy) and intensity rings to hdf5 for further review

    Each fitted ring ("good fit") and failed ring ("no fit") also holds its selection ``status`` and ``reasons``.

    Parameters
    ----------
    iq2d:  ~drtsans.dataobjects.Azimuthal
    rings: ~list
        List of I(Qx, Qy) in rings
    selection: RingSelection
        Fit results and selection outcome of every Q-ring
    output_dir: str
        full path of the output directory

    Returns
    -------

    """
    # annular binning on the full range
    azimuthal_offset = 0.5 * azimuthal_delta
    azimuthal_binning = BinningParams(
        0.0 - azimuthal_offset,
        360.0 - azimuthal_offset,
        bins=int(360.0 / azimuthal_delta),
    )
    q1d = np.sqrt(iq2d.qx**2 + iq2d.qy**2)
    logger.notice(f"[DEBUG] I(Qx, Qy) Q range: {q1d.min()}, {q1d.max()}")
    full_range_annular = bin_annular_into_q1d(iq2d, azimuthal_binning, q1d.min(), q1d.max(), BinningMethod.NOWEIGHT)

    # open
    debug_h5 = h5py.File(os.path.join(output_dir, "auto_wedge_fit.h5"), "w")
    # define h5py string type
    h5_string_type = h5py.special_dtype(vlen=bytes)

    # create 3 groups for result
    ring_group = debug_h5.create_group("rings")
    fit_group = debug_h5.create_group("good fit")
    nofit_group = debug_h5.create_group("no fit")

    # write full range
    group = ring_group.create_group("full range")
    group.create_dataset("phi", data=full_range_annular.phi)
    group.create_dataset("intensity", data=full_range_annular.intensity)

    # write out for each ring
    func_param_dict = dict()
    peak_fit_dict = selection.legacy_dicts()
    ring_fits = {ring_fit.index: ring_fit for ring_fit in selection.ring_fits}

    for index, ring in enumerate(rings):
        # add data to hdf file
        group = ring_group.create_group(f"ring {index}")
        group.create_dataset("phi", data=ring.phi)
        group.create_dataset("intensity", data=ring.intensity)

        # add fitting related information to hdf file
        if peak_fit_dict[index]["error"] is None:
            # good fit
            group = fit_group.create_group(f"ring {index}")
            for param_name in peak_fit_dict[index]:
                # ignore non function parameter
                if param_name in ["fit_function", "error", "used"]:
                    continue

                # create sub dictionary if it does not exist
                if param_name not in func_param_dict:
                    func_param_dict[param_name] = list()

                # add value to dictionary from fit result dictionary
                param_value, param_error = peak_fit_dict[index][param_name]
                func_param_dict[param_name].append([index, param_value, param_error])
        else:
            # no fit
            group = nofit_group.create_group(f"ring {index}")
            # error message
            string_data_set = group.create_dataset("error", (1,), dtype=h5_string_type)
            string_data_set[0] = peak_fit_dict[index]["error"]

        # outcome of the ring selection
        for name, text in [
            ("status", ring_fits[index].status.value),
            ("reasons", "\n".join(ring_fits[index].reasons)),
        ]:
            string_data_set = group.create_dataset(name, (1,), dtype=h5_string_type)
            string_data_set[0] = text

        # fit function
        if "fit_function" in peak_fit_dict[index]:
            function_data_set = group.create_dataset("function", (1,), dtype=h5_string_type)
            function_data_set[0] = peak_fit_dict[index]["fit_function"]

    # add peak fitting result
    for param_name, param_value in func_param_dict.items():
        # form data set
        data_set = np.array(param_value)
        fit_group.create_dataset(param_name, data=data_set)

    # close
    debug_h5.close()


def getWedgeSelection(
    data2d,
    q_min,
    q_delta,
    q_max,
    azimuthal_delta,
    peak_width=0.25,
    background_width=1.5,
    signal_to_noise_min=10.0,
    peak_search_window_size_factor=0.6,
    debug_dir="/tmp/",
    auto_wedge_phi_min=0,
    auto_wedge_phi_max=360,
    auto_symmetric_wedges=False,
    contrast_min=0.2,
    center_tolerance=0.5,
) -> List[List[Tuple[float, float]]]:
    r"""
    Calculate azimuthal binning ranges automatically based on finding peaks in the annular ring. The
    output of this is intended to be used in :py:func:`~drtsans.iq.select_i_of_q_by_wedge`.

    Parameters
    ==========
    data2d: ~drtsans.dataobjects.Azimuthal
        The 2D data object containing the azimuthal data.
    q_min: float
        The left bin boundary for the first Q-bin.
    q_delta: float
        The size of the bins in Q.
    q_max: float
        The left bin boundary for the last Q-bin.
    azimuthal_delta: float
        The size of the bins in azimuthal angle.
    peak_width: float
        Percent of full-width-half-max (FWHM) of the peak to define the signal to be within when
        determining the final range for azimuthal binning.
    background_width: float
        Percent of full-width-half-max (FWHM) of the peak to define the background between peaks
        to be within when determining the final range for azimuthal binning. If the peaks are so wide that the
        background wedges would be empty, it is reduced, with a warning, until each background wedge is as wide
        as the peak wedges.
    signal_to_noise_min: float
        Minimum significance H/sigma_H (fitted peak height over its uncertainty) of each peak of a seed ring.
    peak_search_window_size_factor: float
        Factor of 360 / (num peaks) to construct the search range for wedge peak.
    debug_dir: str
        Full path of the output directory for debugging output files
    auto_wedge_phi_min: float
        minimum angle to find a wedge for subsequent symmetrisation. Used if symmetric_wedges is true
    auto_wedge_phi_max: float
        maximum angle to find a wedge for subsequent symmetrisation. Used if symmetric_wedges is true
    auto_symmetric_wedges: bool
        find a wedge between auto_wedge_phi_min and auto_wedge_phi_max, then symmetrize the findings to the remaining
        phi angle domain
    contrast_min: float
        Minimum contrast H/A0 (fitted peak height over the flat background) of each peak of a seed ring.
    center_tolerance: float
        Maximum deviation of a peak center from the consensus orientation of the seed rings, as a fraction
        of the FWHM of the peak, in (0, 1].

    Returns
    =======
    list
        A two-item list. The first item contains the azimuthal ranges for the two wedges
        that enclose the intensity maxima. The second item contains the azimuthal ranges for the
        background regions associated to the wedges.
        The first item is `[(wedge1\_phimin, wedge1\_phimax), (wedge2\_phimin, wedge2\_phimax)]`,
        and the second item is `[(back1\_phimin, back1\_phimax), (back2\_phimin, back2\_phimax)]]`.

    Raises
    ======
    AutoWedgeError
        If no orientation can be determined from the Q-rings. The diagnostic files are written beforehand.
    """
    criteria = AcceptanceCriteria.from_settings(
        signal_to_noise_min, contrast_min, center_tolerance, azimuthal_delta, maxchisq=1000.0
    )
    if signal_to_noise_min < RESCUE_SIGNIFICANCE_MIN:
        logger.warning(
            f"autoWedgeSignalToNoiseMin = {signal_to_noise_min:g} is the minimum significance H/sigma_H of each "
            "fitted peak, no longer the ratio mean/min of I(phi). Values below "
            f"{RESCUE_SIGNIFICANCE_MIN:g} accept weak peaks and lower the rescue threshold; the default is 10."
        )
    # Bin azimuthal
    q, azimuthal_rings = _binInQAndAzimuthal(
        data2d,
        q_min=q_min,
        q_max=q_max,
        q_delta=q_delta,
        azimuthal_delta=azimuthal_delta,
    )
    selection = _fitQAndAzimuthal(
        azimuthal_rings,
        q_bins=q,
        criteria=criteria,
        azimuthal_start=110.0,
        peak_search_window_size_factor=peak_search_window_size_factor,
        auto_wedge_phi_min=auto_wedge_phi_min,
        auto_wedge_phi_max=auto_wedge_phi_max,
        auto_symmetric_wedges=auto_symmetric_wedges,
    )
    # Export fitting result before any failure is raised, so that the user can inspect the fits
    _export_to_h5(
        iq2d=data2d,
        rings=azimuthal_rings,
        azimuthal_delta=azimuthal_delta,
        selection=selection,
        output_dir=debug_dir,
    )
    _plot_fit_results(azimuthal_rings, selection, debug_dir, auto_symmetric_wedges)

    if not selection.succeeded:
        message = _failure_message(selection, debug_dir)
        logger.error(message)
        raise AutoWedgeError(message, selection)
    logger.notice(f"Automatic wedge detection, Q-ring selection:\n{_format_ring_table(selection)}")
    # the acceptance criteria limit the FWHM of every accepted peak to FWHM_MAX
    center_vec, fwhm_vec = selection.centers, selection.fwhms
    logger.notice(f"Wedge peak centers (deg): {center_vec}, FWHMs (deg): {fwhm_vec}")

    # wide peaks may leave no room for the background wedges at the requested width; reduce it then
    effective_width = _effective_background_width(center_vec, fwhm_vec, background_width, peak_width)
    if effective_width != background_width and effective_width <= peak_width:
        message = (
            f"Automatic wedge detection failed: peaks of FWHM {fwhm_vec[0]:.1f} / {fwhm_vec[1]:.1f} deg leave no "
            f"room for background wedges outside the peak wedges (autoWedgePeakWidth = {peak_width:g}).\n"
            "  Restrict autoWedgeQmin / autoWedgeQmax to the Q-rings with narrower peaks, or set "
            "WedgeMinAngles/WedgeMaxAngles manually.\n\n" + _format_ring_table(selection)
        )
        logger.error(message)
        raise AutoWedgeError(message, selection)
    if background_width <= peak_width:
        logger.warning(
            f"autoWedgeBackgroundWidth = {background_width:g} is not above autoWedgePeakWidth = {peak_width:g}: "
            "the background wedges overlap the peak wedges."
        )
    elif effective_width != background_width:
        logger.warning(
            f"autoWedgeBackgroundWidth reduced from {background_width:g} to {effective_width:.3g}: with peaks of FWHM "
            f"{fwhm_vec[0]:.1f} / {fwhm_vec[1]:.1f} deg, the background wedges would be narrower than the peak "
            "wedges. Each background wedge now has the angular width of the peak wedges."
        )
        background_width = effective_width

    # convert to min and max ranges
    min_vec, max_vec = [], []

    # azimuthal range for the first wedge enclosing an intensity maximum
    min_vec.append(center_vec[0] - peak_width * fwhm_vec[0])
    max_vec.append(center_vec[0] + peak_width * fwhm_vec[0])

    # azimuthal range enclosing the background for the first wedge
    min_vec.append(center_vec[0] + background_width * fwhm_vec[0])
    max_vec.append(center_vec[1] - background_width * fwhm_vec[1])

    # azimuthal range for the second wedge enclosing an intensity maximum
    min_vec.append(center_vec[1] - peak_width * fwhm_vec[1])
    max_vec.append(center_vec[1] + peak_width * fwhm_vec[1])

    # azimuthal range enclosing the background for the second wedge
    min_vec.append(center_vec[1] + background_width * fwhm_vec[1])
    max_vec.append(center_vec[0] - background_width * fwhm_vec[0])

    # clean up the data to be in the form expected by select_i_of_q_by_wedge
    min_vec = np.array(min_vec)
    max_vec = np.array(max_vec)

    min_vec[min_vec < -90.0] += 360.0
    max_vec[max_vec < -90.0] += 360.0

    min_vec[min_vec > 270.0] -= 360.0
    max_vec[max_vec > 270.0] -= 360.0

    # rearrange the ranges: first the wedges, then the backgrounds
    raw_wedges = list(zip(min_vec, max_vec))
    summing_wedges = []
    for i in range(len(raw_wedges) // 2):  # iterate over half the list
        summing_wedges.append([raw_wedges[i], raw_wedges[i + 2]])

    return summing_wedges


def _binInQAndAzimuthal(data, q_min, q_delta, q_max, azimuthal_delta):
    """This function bins the data in Qmod and azimuthal accoring to the supplied parameters. The maximum
    azimuthal is 540deg to allow for finding a peak at/near azimuthal=0deg.

    Parameters
    ==========
    data: ~drtsans.dataobjects.IQazimuthal
    q_min: float
        The left bin boundary for the first Q-bin
    q_delta: float
        The size of the bins in Q
    q_max: float
        The left bin  boundary for the last Q-bin
    azimuthal_delta: float
        The size of the bins in azimuthal

    Results
    =======
    tuple
        Histogram of ```(intensity, error, azimuthal_bins, q_bins)```
    """
    # Export information for Q
    data_q_vec = np.sqrt(data.qx**2 + data.qy**2)
    logger.notice(f"Raw I(Q). Q range: {data_q_vec.min()}, {data_q_vec.max()}")
    # the bonus two steps is to get the end-point in the array
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)
    logger.notice(f"1D Q bins: {q_bins}")

    # create azimuthal binning BinningParams takes number of steps
    azimuthal_offset = 0.5 * azimuthal_delta
    azimuthal_binning = BinningParams(
        0.0 - azimuthal_offset,
        360.0 - azimuthal_offset,
        bins=int(360.0 / azimuthal_delta),
    )
    # create the I(azimuthal) for each q-ring
    data_of_q_rings = []
    # debugging output file
    for qmin_ring, qmax_ring in zip(q_bins[:-1], q_bins[1:]):
        # bin into I(azimuthal)
        i1d_annular = bin_annular_into_q1d(data, azimuthal_binning, qmin_ring, qmax_ring, BinningMethod.NOWEIGHT)

        # Create a copy of the arrays with the 360->540deg region repeated
        # ignore - delta_mod_q wavelength
        phi_new = determine_1d_linear_bins(
            x_min=0.0,
            x_max=540.0 + azimuthal_delta,
            bins=1 + int(540.0 / azimuthal_delta),
        ).centers
        num_orig_bins = i1d_annular.phi.size
        num_repeated_bins = phi_new.size - num_orig_bins

        intensity_new = np.zeros(phi_new.size)
        intensity_new[:num_orig_bins] = i1d_annular.intensity
        intensity_new[-1 * num_repeated_bins :] = i1d_annular.intensity[:num_repeated_bins]

        error_new = np.zeros(phi_new.size)
        error_new[:num_orig_bins] = i1d_annular.error
        error_new[-1 * num_repeated_bins :] = i1d_annular.error[:num_repeated_bins]

        i1d_annular = I1DAnnular(intensity=intensity_new, error=error_new, phi=phi_new)

        # append to the list of spectra
        data_of_q_rings.append(i1d_annular)

    return q_bins, data_of_q_rings


def _estimatePeakParameters(intensity, azimuthal, azimuthal_start, window_half_width):
    """Estimate the peak parameters by determining a window around a bright point in the data then using
    moment calculations to estimate the parameters for a Gaussian that approximates the actual peak.
    This is done to aid the fitting which does better with better starting values.

    Parameters
    ==========
    intensity: numpy.ndarray
        Array of intensities. This must not have nans in it.
    azimuthal: numpy.ndarray
        Array of azimuthal angles centers
    azimuthal_start: float
        Starting guess of peak center
    window_half_width: float
        The window used is this amount on either side of what is determined to be the peak center

    Results
    =======
    tuple
        ``(intensity, center, sigma)`` where intensity is the full height including background
    """
    # Look for the highest point in a section of the data. This is an iterative approach that starts with a window
    # centered at `azimuthal_start`, the repeats until the maximum within the window doesn't move more than 1deg.
    azimuthal_new = azimuthal_start  # where to search around
    # azimuthal_last = azimuthal_start  # last known value
    while True:
        # determine new windows staying at least 90.deg inside the edges
        window_min = np.max((azimuthal_new - window_half_width, azimuthal.min() + 90.0))
        window_max = np.min((azimuthal_new + window_half_width, azimuthal.max() - 90.0))

        # create a search window around azimuthal_new
        left_index = azimuthal.searchsorted(window_min, side="right")
        right_index = azimuthal.searchsorted(window_max, side="right")
        if right_index <= left_index:
            raise RuntimeError(f"no data in the peak search window around {azimuthal_new:.1f} deg")
        # the highest value in the window
        max_value = intensity[left_index:right_index].max()
        # where that is in the window
        max_index = np.where(intensity[left_index:right_index] == max_value)[0].max() + left_index
        # update values
        azimuthal_last = azimuthal_new
        azimuthal_new = azimuthal[max_index]

        # stop searching if the value hasn't changed by less than one degree
        if abs(azimuthal_new - azimuthal_last) < 1.0:
            break
    # output
    logger.debug(
        f"[WEDGE FIT] azimuthal: {azimuthal_new}, {azimuthal_last} with left and right as {left_index}, {right_index}"
    )

    # now use the first two moments of the data within the window to give an improved center position (first moment)
    # and width (derived from second moment)

    # the position of the center of the peak is the first moment of the data. "mean" can be thought of as the
    # center of mass of the peak in azimuthal angle.
    mean = np.sum(intensity[left_index:right_index] * azimuthal[left_index:right_index]) / np.sum(
        intensity[left_index:right_index]
    )

    # the fit uses sigma rather than fwhm
    # calculate the second moment about the mean as an approximation to a Gaussian's "sigma" parameter
    sigma = np.sum(intensity[left_index:right_index] * np.square(azimuthal[left_index:right_index] - mean)) / np.sum(
        intensity[left_index:right_index]
    )
    sigma = np.sqrt(sigma)

    return max_value, mean, sigma


def _fitSpectrum(
    spectrum,
    q_value,
    azimuthal_start,
    peak_search_window_size_factor,
    verbose=True,
    auto_wedge_phi_min=0,
    auto_wedge_phi_max=360,
    auto_symmetric_wedges=False,
):
    """Extract the peak fit parameters for the data. This is done by observing where 2 maxima are in the
    spectrum then fitting for the peak parameters. This makes the assumption that the two peaks are 180deg
    apart.

    Parameters
    ----------
    spectrum:
    q_value: float
        center (of q_min and q_max) for given annular binned Q-ring
    azimuthal_start: float
        First position to look for peaks around
    peak_search_window_size_factor: float
        Factor of 360 / (num peaks) to construct the search range for wedge peak
    verbose: bool
        Flag to output fitting information
    auto_wedge_phi_min: float
        minimum angle to find a wedge for subsequent symmetrisation. Used if symmetric_wedges is true
    auto_wedge_phi_max: float
        maximum angle to find a wedge for subsequent symmetrisation. Used if symmetric_wedges is true
    auto_symmetric_wedges: bool
        find a wedge between auto_wedge_phi_min and auto_wedge_phi_max, then symmetrize the findings to the remaining
        azimuthal angle domain

    Returns
    -------
    dict
        dict[name] = (value, error) where all of the fit parameters are converted.
        f0 is background, then f1...fn are the fitted peaks

    """
    # define a default window size based on the number of peaks the function supports
    # currently only two peaks that are approximately 180deg apart is supported
    # window_factor = 0.6  # default is 0.6 about 108 degree with 2 peaks .. for strong anisotropic: 0.1
    peak_search_window_size = peak_search_window_size_factor * 180
    if auto_symmetric_wedges:
        if auto_wedge_phi_min < 0.0:
            auto_wedge_phi_min += 360.0
            auto_wedge_phi_max += 360.0

        azimuthal_start = (auto_wedge_phi_min + auto_wedge_phi_max) / 2.0
        peak_search_window_size = (auto_wedge_phi_max - auto_wedge_phi_min) / 2.0
    logger.information(
        f"[WEDGE] Fixed window size = {peak_search_window_size}"
        f"from factor {peak_search_window_size_factor} Number of peaks = 2"
    )

    # filter out the nans
    mask = np.logical_not(np.isnan(spectrum.intensity))
    if np.sum(mask) < 10:  # do not allow fitting less points than there are parameters
        raise RuntimeError("Less than 8 points being fit with 7 parameters (found {} points)".format(np.sum(mask)))

    # first estimate background as minimum value
    # this will be subtracted off from found intensities during estimation.
    # Whether the ring is anisotropic is decided after the fit, by the acceptance criteria of the ring selection.
    background = spectrum.intensity[mask].min()

    # start of what will eventually be the fit function by specifying the background
    function = ["name=FlatBackground,A0={}".format(background)]

    # template for describing initial peak guess
    # Sigma is constrained positive so that the reported FWHM (proportional to Sigma) is never negative
    gaussian_str = "name=Gaussian,Height={},PeakCentre={},Sigma={},constraints=(0<Height,0<Sigma)"

    # guess where one peak might be, start with a window of WINDOW_SIZE each side around 110
    intensity_peak, azimuthal_first, sigma = _estimatePeakParameters(
        spectrum.intensity[mask],
        spectrum.phi[mask],
        azimuthal_start=azimuthal_start,
        window_half_width=peak_search_window_size,
    )
    function.append(gaussian_str.format(intensity_peak - background, azimuthal_first, sigma))

    if not auto_symmetric_wedges:
        # assume the other peak is 180 degrees away
        azimuthal_start = azimuthal_first + 180.0
        intensity_peak, azimuthal_second, sigma = _estimatePeakParameters(
            spectrum.intensity[mask],
            spectrum.phi[mask],
            azimuthal_start=azimuthal_start,
            window_half_width=peak_search_window_size,
        )
        function.append(gaussian_str.format(intensity_peak - background, azimuthal_second, sigma))

    # create workspace version of data
    # this includes the nans so `Fit` has to be told to ignore them
    q_azimuthal_workspace = spectrum.to_workspace()

    # fit the positions of the two suspected peaks
    fit_workspace_prefix = mtd.unique_hidden_name()
    fit_function = ";".join(function)
    try:
        if auto_symmetric_wedges:
            fitresult = Fit(
                Function=fit_function,
                InputWorkspace=q_azimuthal_workspace,
                Output=fit_workspace_prefix,
                StartX=spectrum.phi.min() + auto_wedge_phi_min,
                EndX=spectrum.phi.min() + auto_wedge_phi_max,
                OutputParametersOnly=True,
                IgnoreInvalidData=True,
            )
        else:
            fitresult = Fit(
                Function=fit_function,
                InputWorkspace=q_azimuthal_workspace,
                Output=fit_workspace_prefix,
                StartX=spectrum.phi.min() + 90.0,
                EndX=spectrum.phi.min() + 90.0 + 360.0,
                OutputParametersOnly=True,
                IgnoreInvalidData=True,
            )
    except RuntimeError as e:
        raise RuntimeError(f"Failed to fit Q={q_value:.4g} 1/A with fit function {fit_function}") from e
    finally:
        DeleteWorkspace(q_azimuthal_workspace)

    if fitresult.OutputStatus != "success":
        raise RuntimeError(f"Failed to fit Q={q_value:.4g} 1/A")

    # convert the table into a dict[name] = (value, error)
    result = {"fit_function": fit_function}
    for i in range(fitresult.OutputParameters.rowCount()):
        row = fitresult.OutputParameters.row(i)
        logger.debug(f"[DEBUG-TEST] row: {row} of type {type(row)}")
        name = row["Name"]
        if name.startswith("Cost function"):
            name = "chisq"
        result[name] = (row["Value"], row["Error"])

    # delete fit results
    for label in ["_Parameters", "_NormalisedCovarianceMatrix"]:
        DeleteWorkspace(Workspace=fit_workspace_prefix + label)

    # symmetrize the result
    if auto_symmetric_wedges:
        result["f2.Height"] = result["f1.Height"]
        result["f2.PeakCentre"] = (result["f1.PeakCentre"][0] + 180.0, result["f1.PeakCentre"][1])
        result["f2.Sigma"] = result["f1.Sigma"]

    if verbose:
        logger.notice(f"Fit result: {result}")

    return result


def _fitQAndAzimuthal(
    azimuthal_rings,
    q_bins,
    criteria: AcceptanceCriteria,
    azimuthal_start,
    peak_search_window_size_factor,
    verbose=True,
    auto_wedge_phi_min=0,
    auto_wedge_phi_max=360,
    auto_symmetric_wedges=False,
) -> RingSelection:
    """Fit the azimuthal profile of every Q-ring, then select the rings that determine the wedge angles.

    Parameters
    ==========
    azimuthal_rings: list[~drtsans.dataobjects.I1DAnnular]
        Azimuthal profile I(phi) of each Q-ring
    q_bins: numpy.ndarray
        Array of Q bin boundaries, in 1/Angstrom
    criteria: AcceptanceCriteria
        Thresholds of the ring selection
    azimuthal_start: float
        First position to look for peaks around, in degrees
    peak_search_window_size_factor: float
        Factor of 360 / (num peaks) to construct the search range for wedge peak
    verbose: bool
        Flag to turn on fitting information output
    auto_wedge_phi_min, auto_wedge_phi_max, auto_symmetric_wedges:
        See :py:func:`getWedgeSelection`

    Results
    =======
    RingSelection
        Status and reasons of every Q-ring and, on success, the weighted average peak centers and FWHMs.
        A failed selection is reported through ``RingSelection.failure``, never raised.
    """
    if len(azimuthal_rings) != len(q_bins) - 1:
        raise RuntimeError("must supply q-bin boundaries")

    # change to centers in Q for messages
    q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])

    ring_fits = []
    for index, (spectrum, q_center) in enumerate(zip(azimuthal_rings, q_centers)):
        try:
            fitresult = _fitSpectrum(
                spectrum,
                q_center,
                azimuthal_start=azimuthal_start,
                peak_search_window_size_factor=peak_search_window_size_factor,
                verbose=verbose,
                auto_wedge_phi_min=auto_wedge_phi_min,
                auto_wedge_phi_max=auto_wedge_phi_max,
                auto_symmetric_wedges=auto_symmetric_wedges,
            )
        except (RuntimeError, ValueError) as e:
            ring_fits.append(RingFit.failed(index, float(q_center), str(e)))
        else:
            ring_fits.append(RingFit.from_fit_result(index, float(q_center), fitresult, auto_symmetric_wedges))

    return _select_rings(ring_fits, criteria)
