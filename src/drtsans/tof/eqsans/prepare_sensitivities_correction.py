"""EQSANS-specific sensitivity preparation options."""

from __future__ import annotations

import json
import os
from typing import Optional, Sequence, Union

from mantid.api import IEventWorkspace
from mantid.simpleapi import ConvertToMatrixWorkspace, Integration, NormaliseByCurrent, mtd

from drtsans import __version__
from drtsans import subtract_background
from drtsans.prepare_sensivities_correction import (
    PREPARE_DATA,
    PrepareSensitivityCorrection as PrepareBase,
)
from drtsans.samplelogs import SampleLogs
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

    The sensitivity is computed with the detector-patch method, which uses a
    single flood image. Exactly one flood run is therefore accepted, and each
    direct-beam, dark-current, background, and transmission input must also
    contain exactly one run.

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
    def _json_log_value(value) -> str:
        """Serialize a provenance value for a Mantid string log.

        Parameters
        ----------
        value : object
            JSON-compatible value, or a run/configuration value that can be
            represented as text.

        Returns
        -------
        str
            JSON representation, or an empty string for ``None``.
        """
        if value is None:
            return ""
        return json.dumps(value, default=str, sort_keys=True)

    @staticmethod
    def _single_or_json(values):
        """Represent one value directly and multiple values as JSON.

        Parameters
        ----------
        values : sequence
            Values associated with the configured flood runs.

        Returns
        -------
        object
            The only value in ``values``, or a JSON string containing all
            values.
        """
        if not values:
            return ""
        if len(values) == 1:
            return values[0]
        return PrepareSensitivityCorrection._json_log_value(values)

    @staticmethod
    def _workspace_log_value(workspace, name):
        """Read an optional scalar log from a prepared workspace.

        Parameters
        ----------
        workspace : ~mantid.api.MatrixWorkspace or None
            Workspace containing the reduction logs.
        name : str
            Log name to retrieve.

        Returns
        -------
        object
            Log value, or an empty string when the workspace or log is absent.
        """
        if workspace is None:
            return ""
        logs = SampleLogs(workspace)
        if name not in logs:
            return ""
        value = logs[name].value
        if hasattr(value, "item"):
            value = value.item()
        return value

    def _add_sensitivity_provenance(self, sensitivity_ws, provenance):
        """Record EQSANS preparation settings in the sensitivity workspace.

        Parameters
        ----------
        sensitivity_ws : ~mantid.api.MatrixWorkspace
            Export workspace whose sample logs will contain the provenance.
        provenance : dict
            Values collected by :meth:`execute`: ``beam_centers``,
            ``sensitivity_options``, ``wavelength_min``, and
            ``wavelength_max``.

        Notes
        -----
        Values that can contain several runs or beam centers are stored as JSON
        strings so the exported NeXus file remains self-describing.
        """
        options = provenance["sensitivity_options"]
        beam_centers = provenance["beam_centers"]
        beam_center_x = [float(center[0]) for center in beam_centers]
        beam_center_y = [float(center[1]) for center in beam_centers]
        wavelength_min = provenance["wavelength_min"]
        wavelength_max = provenance["wavelength_max"]

        logs = SampleLogs(sensitivity_ws)

        def insert(name, value, unit=None):
            if value is None:
                value = ""
            logs.insert(name, value, unit=unit)

        insert("drtsans_version", __version__)
        insert(
            "sensitivity_method",
            "moving_detectors" if options["use_moving_detector_method"] else "patch",
        )
        insert("flood_runs", self._json_log_value(self._flood_runs))
        insert("beam_center_runs", self._json_log_value(self._direct_beam_center_runs))
        insert("beam_center_x", self._single_or_json(beam_center_x), "m")
        insert("beam_center_y", self._single_or_json(beam_center_y), "m")
        insert("beam_center_mask_radius", self._beam_center_radius, "mm")
        insert("sample_offset", self._sample_offset, "mm")
        insert("detector_offset", self._detector_offset, "mm")
        insert("scale_components", self._json_log_value(self._scale_components))
        insert("solid_angle_correction", int(bool(self._solid_angle_correction)))
        insert("flux_method", "proton charge")
        insert("dark_current_runs", self._json_log_value(self._dark_current_runs))
        insert("blocked_beam_runs", self._json_log_value(self._blocked_beam_runs))
        insert("background_runs", self._json_log_value(self._empty_cell_runs))
        insert("flood_transmission_runs", self._json_log_value(self._transmission_flood_runs))
        insert("background_transmission_runs", self._json_log_value(self._transmission_empty_cell_runs))
        insert("empty_beam_transmission_run", self._json_log_value(self._transmission_reference_runs))
        insert("theta_dependent_transmission", int(self._transmission_reference_runs is not None))
        insert("low_tof_clip", self._low_tof_clip, "microsecond")
        insert("high_tof_clip", self._high_tof_clip, "microsecond")
        insert("wavelength_min", wavelength_min, "Angstrom")
        insert("wavelength_max", wavelength_max, "Angstrom")
        insert("wavelength_weighting", "counts")
        insert("min_threshold", options["min_threshold"])
        insert("max_threshold", options["max_threshold"])
        insert("poly_order", options["polynomial_order"])
        insert("min_detectors_per_tube", options["min_detectors_per_tube"])
        insert("masked_pixels", self._json_log_value(self._extra_mask_dict))
        insert("default_mask", "" if self._default_mask is None else str(self._default_mask))

    def _capture_sensitivity_input_provenance(self, flood_workspace):
        """Capture wavelength logs before sensitivity export reloads the flood run.

        Parameters
        ----------
        flood_workspace : ~mantid.api.MatrixWorkspace
            Prepared EQSANS flood workspace.

        Returns
        -------
        dict
            The flood's ``wavelength_min`` and ``wavelength_max`` logs, in Angstrom.
        """
        return {
            "wavelength_min": self._workspace_log_value(flood_workspace, "wavelength_min"),
            "wavelength_max": self._workspace_log_value(flood_workspace, "wavelength_max"),
        }

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

    @classmethod
    def _as_single_run_list(cls, runs: Optional[RunInput], name: str) -> Optional[list[Run]]:
        """Normalize a run input that must contain exactly one run.

        The EQSANS sensitivity is computed with the detector-patch method,
        which operates on a single flood image. Every run input therefore
        refers to that one flood measurement.

        Parameters
        ----------
        runs : int, str, sequence of int or str, or None
            Run number, run path, or one-element sequence of runs.
        name : str
            Name of the input, used in validation error messages.

        Returns
        -------
        list of int or str or None
            One-element run list, or ``None`` when no runs were supplied.

        Raises
        ------
        TypeError
            If the input contains an invalid run value.
        ValueError
            If the input does not contain exactly one run.
        """
        normalized = cls._as_run_list(runs, name)
        if normalized is not None and len(normalized) != 1:
            raise ValueError(f"{name} must contain exactly one run for EQSANS; got {len(normalized)}")
        return normalized

    def _validate_run_lists(self) -> None:
        """Validate that transmission runs are configured together.

        Raises
        ------
        ValueError
            If transmission runs are only partially supplied.
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

    def set_flood_runs(self, flood_runs: RunInput) -> None:
        """Set the flood run.

        Parameters
        ----------
        flood_runs : int, str, or sequence of int or str
            Flood run number, run path, or one-element sequence of runs.

        Raises
        ------
        TypeError
            If the flood run has an invalid type.
        ValueError
            If more than one flood run is supplied.
        """
        super().set_flood_runs(self._as_single_run_list(flood_runs, "flood_runs"))

    def set_direct_beam_runs(self, direct_beam_runs: RunInput) -> None:
        """Set the direct-beam run used to find the beam center.

        Parameters
        ----------
        direct_beam_runs : int, str, or sequence of int or str
            Direct-beam run number, run path, or one-element sequence of runs.

        Raises
        ------
        TypeError
            If the direct-beam run has an invalid type.
        ValueError
            If more than one direct-beam run is supplied.
        """
        super().set_direct_beam_runs(self._as_single_run_list(direct_beam_runs, "direct_beam_runs"))

    def set_dark_current_runs(self, dark_current_runs: Optional[RunInput]) -> None:
        """Set the optional dark-current run.

        Parameters
        ----------
        dark_current_runs : int, str, sequence of int or str, or None
            Dark-current run number, run path, one-element sequence of runs,
            or ``None``.

        Raises
        ------
        TypeError
            If the dark-current run has an invalid type.
        ValueError
            If more than one dark-current run is supplied.
        """
        super().set_dark_current_runs(self._as_single_run_list(dark_current_runs, "dark_current_runs"))

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
        """Set the time-of-flight clipping widths.

        The two widths are independent: neither is an absolute TOF limit, and
        ``low_tof_clip`` may exceed ``high_tof_clip``.

        Parameters
        ----------
        low_tof_clip : float
            Width trimmed from the fast-neutron end of each transmitted band,
            after its minimum TOF, in microseconds.
        high_tof_clip : float
            Width trimmed from the slow-neutron end of each transmitted band,
            before its maximum TOF, in microseconds.

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
            Blocked-beam run number, run path, one-element sequence of runs, or
            ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If the input does not contain exactly one run.
        """
        self._blocked_beam_runs = self._as_single_run_list(runs, "blocked_beam_runs")

    def set_empty_cell_runs(self, runs: Optional[RunInput]) -> None:
        """Set optional empty-cell runs.

        Parameters
        ----------
        runs : int, str, sequence of int or str, or None
            Empty-cell run number, run path, one-element sequence of runs, or
            ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If the input does not contain exactly one run.
        """
        self._empty_cell_runs = self._as_single_run_list(runs, "empty_cell_runs")

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
            Flood transmission run number, path, one-element sequence of runs,
            or ``None``.
        transmission_empty_cell_runs : int, str, sequence of int or str, or None
            Empty-cell transmission run number, path, one-element sequence of runs,
            or ``None``.
        transmission_reference_runs : int, str, sequence of int or str, or None
            Reference transmission run number, path, one-element sequence of runs,
            or ``None``.

        Raises
        ------
        TypeError
            If a run has an invalid type.
        ValueError
            If only some transmission runs are supplied, or an input does not
            contain exactly one run.
        """
        self._transmission_flood_runs = self._as_single_run_list(transmission_flood_runs, "transmission_flood_runs")
        self._transmission_empty_cell_runs = self._as_single_run_list(
            transmission_empty_cell_runs, "transmission_empty_cell_runs"
        )
        self._transmission_reference_runs = self._as_single_run_list(
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
    def _resolve_run(run: Run, name: str) -> str:
        """Resolve a configured run value for the EQSANS loader.

        Parameters
        ----------
        run : int or str
            Run number or path to a NeXus file.
        name : str
            Name of the configured run, used in validation errors.

        Returns
        -------
        str
            Instrument-qualified run name or validated run path.

        Raises
        ------
        ValueError
            If the run value is not recognized or a file path does not exist.
        """
        if isinstance(run, str) and run.isdigit():
            run = int(run)
        if isinstance(run, int):
            return f"EQSANS_{run}"
        if isinstance(run, str) and (run.startswith("EQSANS") or os.path.exists(run)):
            return run
        raise ValueError(f"{name} {run!r} is not recognized")

    def _resolve_beam_center_run(self, beam_center_run: Run) -> str:
        """Resolve the direct-beam run with the EQSANS run-name rules.

        Parameters
        ----------
        beam_center_run : int or str
            Run number, instrument-qualified run name, or path to a NeXus file.

        Returns
        -------
        str
            Instrument-qualified run name or validated run path.
        """
        return self._resolve_run(beam_center_run, "direct_beam_run")

    def _resolve_parent_flood_run(self, parent_flood_run: Run) -> str:
        """Resolve the flood run used as the parent of the exported map.

        Parameters
        ----------
        parent_flood_run : int or str
            Run number, instrument-qualified run name, or path to a NeXus file.

        Returns
        -------
        str
            Instrument-qualified run name or validated run path.
        """
        return self._resolve_run(parent_flood_run, "flood_run")

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
            Prepared, proton-charge-normalized histogram workspace with
            wavelength bins.
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
        workspace = NormaliseByCurrent(InputWorkspace=workspace, OutputWorkspace=str(workspace))
        # Without a dark-current run, prepare_data can return an EventWorkspace (when no
        # histogram bin is empty), which background subtraction does not accept. The
        # conversion keeps the wavelength bins and their uncertainties.
        if isinstance(workspace, IEventWorkspace):
            workspace = ConvertToMatrixWorkspace(InputWorkspace=workspace, OutputWorkspace=str(workspace))
        return workspace

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

    def execute(
        self,
        use_moving_detector_method: bool,
        min_threshold: float,
        max_threshold: float,
        output_nexus_name: str,
        enforce_use_nexus_idf: bool = False,
        debug_mode: bool = False,
    ) -> None:
        """Calculate and export the EQSANS sensitivity correction.

        EQSANS sensitivities are calculated with the detector-patch method.
        The moving-detector method requires floods in which the beam-stop
        shadow falls on different pixels, and the common workflow implements
        it only for GPSANS and BIOSANS.

        Parameters
        ----------
        use_moving_detector_method : bool
            Must be ``False``; the moving-detector method is not supported for
            EQSANS.
        min_threshold : float
            Minimum normalized detector response accepted for a good pixel.
        max_threshold : float
            Maximum normalized detector response accepted for a good pixel.
        output_nexus_name : str
            Path to the output processed NeXus file.
        enforce_use_nexus_idf : bool, optional
            Whether to use the instrument definition stored in the NeXus file.
        debug_mode : bool, optional
            Whether to save intermediate flood workspaces for debugging.

        Raises
        ------
        ValueError
            If ``use_moving_detector_method`` is ``True``, or if transmission
            correction is configured without an empty-cell run.
        """
        if use_moving_detector_method:
            raise ValueError("The moving-detector method is not supported for EQSANS; use the detector-patch method")
        # Transmission is applied to the flood and the empty cell before their
        # subtraction, so it requires an empty-cell run. Check this before the
        # beam-center and flood workspaces are loaded.
        if self._transmission_reference_runs is not None and self._empty_cell_runs is None:
            raise ValueError("empty_cell_runs must be configured for transmission correction")
        super().execute(
            use_moving_detector_method,
            min_threshold,
            max_threshold,
            output_nexus_name,
            enforce_use_nexus_idf=enforce_use_nexus_idf,
            debug_mode=debug_mode,
        )
