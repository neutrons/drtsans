"""EQSANS-specific sensitivity preparation options."""

from __future__ import annotations

from typing import Optional, Sequence, Union

from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection as PrepareBase


Run = Union[int, str]
RunInput = Union[Run, Sequence[Run]]


class PrepareSensitivityCorrection(PrepareBase):
    """Prepare EQSANS sensitivity corrections.

    EQSANS-specific preparation options are kept here instead of being added
    to the common sensitivity-preparation interface.  The correction pipeline
    itself will be implemented in later stages.

    Parameters
    ----------
    component : str, optional
        Detector component to use for the sensitivity correction.
    """

    def __init__(self, component: str = "detector1") -> None:
        super().__init__(instrument="EQSANS", component=component)
        self._sample_offset = 0.0
        self._detector_offset = 0.0
        self._low_tof_clip = 500.0
        self._high_tof_clip = 2000.0
        self._blocked_beam_runs: Optional[list[Run]] = None
        self._empty_cell_runs: Optional[list[Run]] = None
        self._transmission_flood_runs: Optional[list[Run]] = None
        self._transmission_empty_cell_runs: Optional[list[Run]] = None
        self._transmission_reference_runs: Optional[list[Run]] = None

    @staticmethod
    def _as_run_list(runs: Optional[RunInput], name: str) -> Optional[list[Run]]:
        """Normalize one run or a sequence of runs to a list.

        Parameters
        ----------
        runs : int, str, sequence of int or str, or None
            Run number, run path, or sequence of runs to normalize.
        name : str
            Name of the input, used in validation error messages.

        Returns
        -------
        list of int or str or None
            Normalized run list, or ``None`` when no runs were supplied.

        Raises
        ------
        TypeError
            If the input is not a run number, path, or sequence of runs, or
            contains an invalid run value.
        ValueError
            If the input sequence is empty.
        """
        if runs is None:
            return None
        if isinstance(runs, (int, str)):
            return [runs]
        try:
            normalized = list(runs)
        except TypeError as exc:
            raise TypeError(f"{name} must be a run number, path, or sequence of runs") from exc
        if not normalized:
            raise ValueError(f"{name} cannot be empty")
        if not all(isinstance(run, (int, str)) for run in normalized):
            raise TypeError(f"{name} must contain only run numbers or paths")
        return normalized

    def _validate_run_lists(self) -> None:
        """Validate and broadcast optional run lists after flood runs are known.

        Raises
        ------
        ValueError
            If transmission runs are only partially supplied, or an auxiliary
            run list does not contain one run or one run per flood run.
        """
        transmission_runs = (
            self._transmission_flood_runs,
            self._transmission_empty_cell_runs,
            self._transmission_reference_runs,
        )
        if any(runs is None for runs in transmission_runs) and any(runs is not None for runs in transmission_runs):
            raise ValueError(
                "transmission_flood_runs, transmission_empty_cell_runs, and "
                "transmission_reference_runs must be supplied together"
            )

        if self._flood_runs is None:
            return

        number_of_flood_runs = len(self._flood_runs)
        attributes = (
            "_blocked_beam_runs",
            "_empty_cell_runs",
            "_transmission_flood_runs",
            "_transmission_empty_cell_runs",
            "_transmission_reference_runs",
        )
        for attribute in attributes:
            runs = getattr(self, attribute)
            if runs is None:
                continue
            if len(runs) == 1 and number_of_flood_runs > 1:
                setattr(self, attribute, runs * number_of_flood_runs)
            elif len(runs) != number_of_flood_runs:
                raise ValueError(
                    f"{attribute[1:]} must contain one run or {number_of_flood_runs} runs; got {len(runs)}"
                )

    def set_flood_runs(self, flood_runs: RunInput) -> None:
        """Set flood runs and validate configured auxiliary runs.

        Parameters
        ----------
        flood_runs : int, str, or sequence of int or str
            Flood run number, run path, or sequence of flood runs.

        Raises
        ------
        TypeError
            If a flood run has an invalid type.
        ValueError
            If configured auxiliary run lists do not match the flood run
            count.
        """
        super().set_flood_runs(flood_runs)
        self._validate_run_lists()

    def set_sample_detector_offsets(self, sample_offset: float, detector_offset: float) -> None:
        """Set sample and detector offsets.

        Parameters
        ----------
        sample_offset : float
            Distance from the sample position, in millimeters.
        detector_offset : float
            Distance from the detector position, in millimeters.
        """
        self._sample_offset = float(sample_offset)
        self._detector_offset = float(detector_offset)

    def set_tof_cuts(self, low_tof_clip: float, high_tof_clip: float) -> None:
        """Set low and high time-of-flight clipping values.

        Parameters
        ----------
        low_tof_clip : float
            Lower time-of-flight limit, in microseconds.
        high_tof_clip : float
            Upper time-of-flight limit, in microseconds.

        Raises
        ------
        ValueError
            If either clipping value is negative.
        """
        low_tof_clip = float(low_tof_clip)
        high_tof_clip = float(high_tof_clip)
        if low_tof_clip < 0.0 or high_tof_clip < 0.0:
            raise ValueError("TOF clipping values must be non-negative")
        self._low_tof_clip = low_tof_clip
        self._high_tof_clip = high_tof_clip

    def set_blocked_beam_runs(self, runs: Optional[RunInput]) -> None:
        """Set optional blocked-beam runs.

        Parameters
        ----------
        runs : int, str, sequence of int or str, or None
            Blocked-beam run number, run path, sequence of runs, or ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If the run list is empty or does not match the flood run count.
        """
        self._blocked_beam_runs = self._as_run_list(runs, "blocked_beam_runs")
        self._validate_run_lists()

    def set_empty_cell_runs(self, runs: Optional[RunInput]) -> None:
        """Set optional empty-cell runs.

        Parameters
        ----------
        runs : int, str, sequence of int or str, or None
            Empty-cell run number, run path, sequence of runs, or ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If the run list is empty or does not match the flood run count.
        """
        self._empty_cell_runs = self._as_run_list(runs, "empty_cell_runs")
        self._validate_run_lists()

    def set_transmission_correction(
        self,
        transmission_flood_runs: Optional[RunInput],
        transmission_empty_cell_runs: Optional[RunInput],
        transmission_reference_runs: Optional[RunInput],
    ) -> None:
        """Set runs used for transmission correction.

        Parameters
        ----------
        transmission_flood_runs : int, str, sequence of int or str, or None
            Flood transmission run number, path, sequence of runs, or
            ``None``.
        transmission_empty_cell_runs : int, str, sequence of int or str, or None
            Empty-cell transmission run number, path, sequence of runs, or
            ``None``.
        transmission_reference_runs : int, str, sequence of int or str, or None
            Reference transmission run number, path, sequence of runs, or
            ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If only some transmission runs are supplied, a run list is empty,
            or a run list does not match the flood run count.
        """
        self._transmission_flood_runs = self._as_run_list(transmission_flood_runs, "transmission_flood_runs")
        self._transmission_empty_cell_runs = self._as_run_list(
            transmission_empty_cell_runs, "transmission_empty_cell_runs"
        )
        self._transmission_reference_runs = self._as_run_list(
            transmission_reference_runs, "transmission_reference_runs"
        )
        self._validate_run_lists()

    def _prepare_data_opts(self, beam_center):
        """Set EQSANS options for flood and auxiliary data loads.

        Parameters
        ----------
        beam_center : list
            Beam-center coordinates passed to the instrument loader.

        Returns
        -------
        dict
            Options for the instrument's ``prepare_data`` function, including
            EQSANS geometry and time-of-flight settings.
        """
        options = super()._prepare_data_opts(beam_center)
        options.update(
            detector_offset=self._detector_offset,
            sample_offset=self._sample_offset,
            low_tof_clip=self._low_tof_clip,
            high_tof_clip=self._high_tof_clip,
        )
        return options

    def _prepare_beam_center_data_opts(self, beam_center):
        """Set options for loading a beam-center workspace.

        Parameters
        ----------
        beam_center : list
            Beam-center coordinates passed to the instrument loader.

        Returns
        -------
        dict
            Base instrument-loader options without EQSANS flood geometry or
            time-of-flight settings.
        """
        return PrepareBase._prepare_data_opts(self, beam_center)
