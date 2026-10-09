"""EQSANS-specific sensitivity preparation options."""

from __future__ import annotations

import os
from typing import Optional, Sequence, Union

from mantid.simpleapi import Integration, NormaliseByCurrent, mtd

from drtsans import subtract_background
from drtsans.prepare_sensivities_correction import (
    PREPARE_DATA,
    PrepareSensitivityCorrection as PrepareBase,
)
from drtsans.tof.eqsans.transmission import calculate_transmission
from drtsans.transmission import apply_transmission_correction


Run = Union[int, str]
RunInput = Union[Run, Sequence[Run]]


class PrepareSensitivityCorrection(PrepareBase):
    """Prepare EQSANS sensitivity corrections.

    EQSANS-specific preparation options are kept here instead of being added
    to the common sensitivity-preparation interface. Flood data is normalized
    by proton charge and corrected for configured blocked-beam and empty-cell
    backgrounds and wavelength-dependent transmission corrections before
    wavelength summation.

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

    @staticmethod
    def _resolve_run(run: Run, name: str) -> Run:
        """Resolve a configured run value for the EQSANS loader.

        Parameters
        ----------
        run : int or str
            Run number or path to a NeXus file.
        name : str
            Name of the configured run, used in validation errors.

        Returns
        -------
        int or str
            Integer run number or validated run path.

        Raises
        ------
        ValueError
            If the run value is not recognized or a file path does not exist.
        """
        if isinstance(run, str) and run.isdigit():
            return int(run)
        if isinstance(run, int):
            return run
        if isinstance(run, str) and (run.startswith("EQSANS") or os.path.exists(run)):
            return run
        raise ValueError(f"{name} {run!r} is not recognized")

    def _prepare_eqsans_workspace(
        self,
        run: Run,
        beam_center,
        dark_current_run,
        name: str,
        solid_angle: Optional[bool] = None,
    ):
        """Load and prepare one EQSANS workspace without wavelength summation.

        Parameters
        ----------
        run : int or str
            Run number or path to a NeXus file.
        beam_center : list
            Beam-center coordinates passed to the EQSANS loader.
        dark_current_run : int, str, or None
            Optional dark-current run associated with the workspace.
        name : str
            Name of the run configuration, used in validation errors.
        solid_angle : bool or None, optional
            Whether to apply solid-angle correction. If ``None``, use the
            preparer's configured value.

        Returns
        -------
        ~mantid.api.MatrixWorkspace
            Prepared, proton-charge-normalized workspace with wavelength bins.
        """
        resolved_run = self._resolve_run(run, name)
        resolved_dark_current = (
            None if dark_current_run is None else self._resolve_run(dark_current_run, "dark_current_run")
        )
        options = dict(
            data=resolved_run,
            scale_components=self.scale_components,
            pixel_calibration=self._apply_calibration,
            mask=self._default_mask,
            btp=self._extra_mask_dict,
            dark_current=resolved_dark_current,
            solid_angle=self._solid_angle_correction if solid_angle is None else solid_angle,
            **self._prepare_data_opts(beam_center),
        )
        if name != "flood_run":
            options["output_workspace"] = mtd.unique_hidden_name()
        workspace = PREPARE_DATA["EQSANS"](**options)
        return NormaliseByCurrent(InputWorkspace=workspace, OutputWorkspace=str(workspace))

    @staticmethod
    def _integrate_workspace(workspace):
        """Sum all wavelength bins in a prepared workspace.

        Parameters
        ----------
        workspace : ~mantid.api.MatrixWorkspace
            Workspace whose wavelength bins should be summed.

        Returns
        -------
        ~mantid.api.MatrixWorkspace
            Workspace containing one integrated bin per detector spectrum.
        """
        if workspace.blocksize() != 1:
            return Integration(InputWorkspace=workspace, OutputWorkspace=str(workspace))
        return workspace

    @staticmethod
    def _delete_workspace(workspace) -> None:
        """Delete a temporary workspace when it is registered.

        Parameters
        ----------
        workspace : ~mantid.api.MatrixWorkspace or None
            Workspace to remove from the Mantid analysis data service.
        """
        if workspace is not None and mtd.doesExist(str(workspace)):
            mtd.remove(str(workspace))

    def _apply_transmission_pair(self, input_workspace, transmission_sample, transmission_reference):
        """Apply a fitted EQSANS transmission calculated from two runs.

        Parameters
        ----------
        input_workspace : ~mantid.api.MatrixWorkspace
            Workspace to correct.
        transmission_sample : ~mantid.api.MatrixWorkspace
            Charge-normalized transmission sample workspace.
        transmission_reference : ~mantid.api.MatrixWorkspace
            Charge-normalized empty-beam transmission workspace.

        Returns
        -------
        ~mantid.api.MatrixWorkspace
            Transmission-corrected input workspace.
        """
        transmission_workspace = calculate_transmission(
            transmission_sample,
            transmission_reference,
        )
        try:
            return apply_transmission_correction(
                input_workspace,
                trans_workspace=transmission_workspace,
                theta_dependent=True,
                output_workspace=str(input_workspace),
            )
        finally:
            self._delete_workspace(transmission_workspace)

    def _flood_run_index(self, flood_run: Run) -> int:
        """Find the configured index for a flood run.

        Parameters
        ----------
        flood_run : int or str
            Flood run currently being prepared.

        Returns
        -------
        int
            Index of ``flood_run`` in the configured flood-run list.
        """
        for index, configured_run in enumerate(self._flood_runs):
            if configured_run == flood_run or str(configured_run) == str(flood_run):
                return index
        raise ValueError(f"Flood run {flood_run!r} is not configured")

    def _prepare_flood_data(self, flood_run, beam_center, dark_current_run):
        """Prepare an EQSANS flood with charge normalization and backgrounds.

        The flood, blocked-beam, and empty-cell workspaces are prepared with
        the same geometry and wavelength cuts. Each is normalized by proton
        charge without using a beam-flux spectrum. Backgrounds are subtracted
        before the flood is integrated over wavelength.

        Parameters
        ----------
        flood_run : int or str
            Flood run number or path to a NeXus file.
        beam_center : list
            Beam-center coordinates used to center all loaded workspaces.
        dark_current_run : int, str, or None
            Optional dark-current run associated with the flood and auxiliary
            workspaces.

        Returns
        -------
        ~mantid.api.MatrixWorkspace
            Corrected flood workspace with one integrated bin per spectrum.
        """
        run_index = self._flood_run_index(flood_run)
        if self._transmission_reference_runs is not None and self._empty_cell_runs is None:
            raise ValueError("empty_cell_runs must be configured for transmission correction")

        flood_workspace = self._prepare_eqsans_workspace(
            flood_run,
            beam_center,
            dark_current_run,
            "flood_run",
        )

        blocked_workspace = None
        empty_cell_workspace = None
        transmission_reference_workspace = None
        transmission_flood_workspace = None
        transmission_empty_cell_workspace = None
        try:
            if self._blocked_beam_runs is not None:
                blocked_workspace = self._prepare_eqsans_workspace(
                    self._blocked_beam_runs[run_index],
                    beam_center,
                    dark_current_run,
                    "blocked_beam_run",
                )
                flood_workspace = subtract_background(
                    flood_workspace,
                    blocked_workspace,
                    output_workspace=str(flood_workspace),
                )

            if self._empty_cell_runs is not None:
                empty_cell_workspace = self._prepare_eqsans_workspace(
                    self._empty_cell_runs[run_index],
                    beam_center,
                    dark_current_run,
                    "empty_cell_run",
                )
                if blocked_workspace is not None:
                    empty_cell_workspace = subtract_background(
                        empty_cell_workspace,
                        blocked_workspace,
                        output_workspace=str(empty_cell_workspace),
                    )

                if self._transmission_reference_runs is not None:
                    transmission_reference_workspace = self._prepare_eqsans_workspace(
                        self._transmission_reference_runs[run_index],
                        beam_center,
                        dark_current_run,
                        "transmission_reference_run",
                        solid_angle=False,
                    )
                    transmission_flood_workspace = self._prepare_eqsans_workspace(
                        self._transmission_flood_runs[run_index],
                        beam_center,
                        dark_current_run,
                        "transmission_flood_run",
                        solid_angle=False,
                    )
                    transmission_empty_cell_workspace = self._prepare_eqsans_workspace(
                        self._transmission_empty_cell_runs[run_index],
                        beam_center,
                        dark_current_run,
                        "transmission_empty_cell_run",
                        solid_angle=False,
                    )
                    flood_workspace = self._apply_transmission_pair(
                        flood_workspace,
                        transmission_flood_workspace,
                        transmission_reference_workspace,
                    )
                    empty_cell_workspace = self._apply_transmission_pair(
                        empty_cell_workspace,
                        transmission_empty_cell_workspace,
                        transmission_reference_workspace,
                    )

                flood_workspace = subtract_background(
                    flood_workspace,
                    empty_cell_workspace,
                    output_workspace=str(flood_workspace),
                )
        finally:
            self._delete_workspace(empty_cell_workspace)
            self._delete_workspace(blocked_workspace)
            self._delete_workspace(transmission_flood_workspace)
            self._delete_workspace(transmission_empty_cell_workspace)
            self._delete_workspace(transmission_reference_workspace)

        return self._integrate_workspace(flood_workspace)

    def _apply_transmission_correction(self, input_workspace, *args):
        """Return the workspace after EQSANS transmission correction.

        EQSANS transmission is applied to the flood and empty-cell workspaces
        in :meth:`_prepare_flood_data`, while they still contain wavelength
        bins. The common executor calls this hook afterward for historical
        compatibility; applying a second correction there would be incorrect.

        Parameters
        ----------
        input_workspace : ~mantid.api.MatrixWorkspace
            Already transmission-corrected workspace.
        *args
            Legacy common-preparer transmission arguments, intentionally
            ignored.

        Returns
        -------
        ~mantid.api.MatrixWorkspace
            The unchanged input workspace.
        """
        return input_workspace

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
