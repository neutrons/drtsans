"""Prepare an EQSANS sensitivity correction map from a TOML configuration."""

import argparse
import copy
import math
import tomllib
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union, cast

from drtsans.tof.eqsans.prepare_sensitivities_correction import PrepareSensitivityCorrection


Run = Union[int, str]
RunInput = Union[Run, Sequence[Run]]
OptionalRunInput = Optional[RunInput]
ScaleComponents = dict[str, list[float]]

DOCUMENTATION_URL = "https://drtsans.readthedocs.io/latest/user/sensitivity_eqsans.html"


# Run inputs contain exactly one run: an integer run number, a string run
# identifier/path, or a one-element sequence of either. Examples: 12345,
# "EQSANS_12345", "/path/to/run.nxs", or (12345,). The EQSANS detector-patch
# method computes the sensitivity from a single flood run. The flood and
# direct-beam runs are required and have no default value.
FLOOD_RUNS: OptionalRunInput = None
DIRECT_BEAM_RUNS: OptionalRunInput = None

# Optional auxiliary run inputs use the same forms as the flood run. Leave
# these as None when the correction is not needed.
DARK_CURRENT_RUNS: OptionalRunInput = None
BLOCKED_BEAM_RUNS: OptionalRunInput = None
EMPTY_CELL_RUNS: OptionalRunInput = None

# Optional wavelength-dependent transmission runs use the same run-input forms.
# Configure all three together: transmission flood, transmission empty-cell,
# and transmission reference.
TRANSMISSION_FLOOD_RUNS: OptionalRunInput = None
TRANSMISSION_EMPTY_CELL_RUNS: OptionalRunInput = None
TRANSMISSION_REFERENCE_RUNS: OptionalRunInput = None

# EQSANS geometry and time-of-flight options.
SAMPLE_OFFSET = 0.0  # mm
DETECTOR_OFFSET = 0.0  # mm
# TOF clipping widths are independent, non-negative amounts in microseconds.
# LOW_TOF_CLIP trims the fast-neutron end of each transmitted band (after its
# minimum TOF); HIGH_TOF_CLIP trims the slow-neutron end (before its maximum).
LOW_TOF_CLIP = 500.0  # microseconds
HIGH_TOF_CLIP = 2000.0  # microseconds

# Masks and detector preparation.
# Universal detector mask.
# Input type: str mask-file path, Mantid MaskWorkspace, list of detector IDs,
# or None. Example: "/path/to/eqsans_mask.xml".
UNIVERSAL_MASK = None
# Additional pixels to mask.
# Input type: str pixel ranges, or None. Example: "1-18,239-256".
MASKED_PIXELS = None
# Radius of the beam-center mask.
# Input type: float in millimeters. Example: 65.0.
MASK_BEAM_CENTER_RADIUS = 65.0
# Apply pixel-width and pixel-height calibration.
# Input type: bool or str calibration-file path. Use False to disable
# calibration, True for the default calibration, or a path for a specific file.
PIXEL_CALIBRATION = False
# Scale detector components along their X, Y, and Z axes.
# Reduction configuration equivalent: configuration.scaleComponents.
# Input type: dict[str, list[float]] or None. Each list contains positive
# X, Y, and Z scale factors. For EQSANS, detector1 is the supported component.
# Example: {"detector1": [1.0, 1.1, 1.0]}.
SCALE_COMPONENTS = None
# Apply solid-angle correction to the flood data.
# Input type: bool. Example: True.
SOLID_ANGLE_CORRECTION = True

# Sensitivity calculation.
# Select the sensitivity-preparation method.
# Input type: bool. Must be False: EQSANS uses the detector-patch method, and
# the moving-detector method is not supported. Example: False.
MOVING_DETECTORS = False
# Minimum accepted normalized detector response.
# Input type: float. Pixels below this value are masked. Example: 0.5.
MIN_THRESHOLD = 0.5
# Maximum accepted normalized detector response.
# Input type: float. Pixels above this value are masked. Example: 2.0.
MAX_THRESHOLD = 2.0


DEFAULT_CONFIGURATION: Mapping[str, Any] = {
    "FLOOD_RUNS": FLOOD_RUNS,
    "DIRECT_BEAM_RUNS": DIRECT_BEAM_RUNS,
    "DARK_CURRENT_RUNS": DARK_CURRENT_RUNS,
    "BLOCKED_BEAM_RUNS": BLOCKED_BEAM_RUNS,
    "EMPTY_CELL_RUNS": EMPTY_CELL_RUNS,
    "TRANSMISSION_FLOOD_RUNS": TRANSMISSION_FLOOD_RUNS,
    "TRANSMISSION_EMPTY_CELL_RUNS": TRANSMISSION_EMPTY_CELL_RUNS,
    "TRANSMISSION_REFERENCE_RUNS": TRANSMISSION_REFERENCE_RUNS,
    "SAMPLE_OFFSET": SAMPLE_OFFSET,
    "DETECTOR_OFFSET": DETECTOR_OFFSET,
    "LOW_TOF_CLIP": LOW_TOF_CLIP,
    "HIGH_TOF_CLIP": HIGH_TOF_CLIP,
    "UNIVERSAL_MASK": UNIVERSAL_MASK,
    "MASKED_PIXELS": MASKED_PIXELS,
    "MASK_BEAM_CENTER_RADIUS": MASK_BEAM_CENTER_RADIUS,
    "PIXEL_CALIBRATION": PIXEL_CALIBRATION,
    "SCALE_COMPONENTS": SCALE_COMPONENTS,
    "SOLID_ANGLE_CORRECTION": SOLID_ANGLE_CORRECTION,
    "MOVING_DETECTORS": MOVING_DETECTORS,
    "MIN_THRESHOLD": MIN_THRESHOLD,
    "MAX_THRESHOLD": MAX_THRESHOLD,
}

_RUN_OPTIONS = {
    "FLOOD_RUNS",
    "DIRECT_BEAM_RUNS",
    "DARK_CURRENT_RUNS",
    "BLOCKED_BEAM_RUNS",
    "EMPTY_CELL_RUNS",
    "TRANSMISSION_FLOOD_RUNS",
    "TRANSMISSION_EMPTY_CELL_RUNS",
    "TRANSMISSION_REFERENCE_RUNS",
}
_REQUIRED_RUN_OPTIONS = ("FLOOD_RUNS", "DIRECT_BEAM_RUNS")
_OPTIONAL_RUN_OPTIONS = _RUN_OPTIONS - set(_REQUIRED_RUN_OPTIONS)
_FLOAT_OPTIONS = {
    "SAMPLE_OFFSET",
    "DETECTOR_OFFSET",
    "LOW_TOF_CLIP",
    "HIGH_TOF_CLIP",
    "MASK_BEAM_CENTER_RADIUS",
    "MIN_THRESHOLD",
    "MAX_THRESHOLD",
}
_BOOL_OPTIONS = {"SOLID_ANGLE_CORRECTION", "MOVING_DETECTORS"}


def _normalize_run(value: Any, option_name: str) -> Run:
    """Normalize one run value from a configuration file."""
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ValueError(f"{option_name} must contain integers or strings")
    if isinstance(value, int):
        return value

    value = value.strip()
    if not value:
        raise ValueError(f"{option_name} cannot contain an empty string")
    try:
        return int(value)
    except ValueError:
        return value


def _normalize_run_input(value: Any, option_name: str, optional: bool = False) -> Optional[tuple[Run, ...]]:
    """Normalize a scalar or sequence of run values to a tuple."""
    if value is None:
        if optional:
            return None
        raise ValueError(f"{option_name} cannot be null")

    values = value if isinstance(value, (list, tuple)) else [value]
    if not values:
        raise ValueError(f"{option_name} cannot be empty")
    return tuple(_normalize_run(item, option_name) for item in values)


def _normalize_float(value: Any, option_name: str) -> float:
    """Normalize a finite numeric configuration value to a float."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{option_name} must be a number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{option_name} must be finite")
    return value


def _normalize_bool(value: Any, option_name: str) -> bool:
    """Validate a configuration boolean."""
    if not isinstance(value, bool):
        raise ValueError(f"{option_name} must be a boolean")
    return value


def _normalize_scale_components(value: Any) -> Optional[ScaleComponents]:
    """Normalize component scale factors from a TOML inline table."""
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("SCALE_COMPONENTS must be a dictionary")

    normalized: ScaleComponents = {}
    for component, factors in value.items():
        if not isinstance(component, str) or not isinstance(factors, list) or len(factors) != 3:
            raise ValueError("SCALE_COMPONENTS values must be three-number lists")
        normalized_factors = [_normalize_float(factor, f"SCALE_COMPONENTS[{component}]") for factor in factors]
        if any(factor <= 0.0 for factor in normalized_factors):
            raise ValueError("SCALE_COMPONENTS factors must be positive")
        normalized[component] = normalized_factors
    return normalized


def _normalize_option(option_name: str, value: Any) -> Any:
    """Normalize one supported configuration option."""
    if option_name in _RUN_OPTIONS:
        return _normalize_run_input(value, option_name, optional=option_name in _OPTIONAL_RUN_OPTIONS)
    if option_name in _FLOAT_OPTIONS:
        return _normalize_float(value, option_name)
    if option_name in _BOOL_OPTIONS:
        return _normalize_bool(value, option_name)
    if option_name == "SCALE_COMPONENTS":
        return _normalize_scale_components(value)
    if option_name == "UNIVERSAL_MASK":
        if value is None or isinstance(value, str):
            return value
        if isinstance(value, list) and all(
            isinstance(detector_id, int) and not isinstance(detector_id, bool) for detector_id in value
        ):
            return value
        raise ValueError("UNIVERSAL_MASK must be a string, a list of detector IDs, or omitted")
    if option_name == "MASKED_PIXELS":
        if value is None or isinstance(value, str):
            return value
        raise ValueError("MASKED_PIXELS must be a string or omitted")
    if option_name == "PIXEL_CALIBRATION":
        if isinstance(value, bool) or isinstance(value, str):
            return value
        raise ValueError("PIXEL_CALIBRATION must be a boolean or a file path")
    raise ValueError(f"Unknown configuration option: {option_name}")


def merge_configuration(overrides: Mapping[str, Any]) -> dict[str, Any]:
    """Merge and validate TOML settings with the script defaults.

    Parameters
    ----------
    overrides
        Configuration values read from the TOML file.

    Returns
    -------
    dict[str, Any]
        A normalized, complete configuration.

    Raises
    ------
    ValueError
        If an option is unknown, a required option is missing, or an option
        has an invalid value.
    """
    unknown_options = set(overrides) - set(DEFAULT_CONFIGURATION)
    if unknown_options:
        names = ", ".join(sorted(unknown_options))
        raise ValueError(f"Unknown configuration option(s): {names}")
    missing_options = [name for name in _REQUIRED_RUN_OPTIONS if name not in overrides]
    if missing_options:
        names = ", ".join(missing_options)
        raise ValueError(f"Missing required configuration option(s): {names}")

    configuration = copy.deepcopy(dict(DEFAULT_CONFIGURATION))
    for option_name, value in overrides.items():
        configuration[option_name] = _normalize_option(option_name, value)

    transmission_values = [
        configuration[name]
        for name in (
            "TRANSMISSION_FLOOD_RUNS",
            "TRANSMISSION_EMPTY_CELL_RUNS",
            "TRANSMISSION_REFERENCE_RUNS",
        )
    ]
    if any(value is not None for value in transmission_values) and not all(
        value is not None for value in transmission_values
    ):
        raise ValueError("TRANSMISSION_* options must be supplied together")
    if transmission_values[0] is not None and configuration["EMPTY_CELL_RUNS"] is None:
        raise ValueError("EMPTY_CELL_RUNS is required when TRANSMISSION_* options are supplied")

    for option_name in sorted(_RUN_OPTIONS):
        runs = configuration[option_name]
        if runs is not None and len(runs) != 1:
            raise ValueError(f"{option_name} must contain exactly one run")

    if configuration["MOVING_DETECTORS"]:
        raise ValueError("MOVING_DETECTORS must be false; the moving-detector method is not supported for EQSANS")
    for option_name in ("LOW_TOF_CLIP", "HIGH_TOF_CLIP"):
        if configuration[option_name] < 0.0:
            raise ValueError(f"{option_name} must be non-negative")
    if configuration["MIN_THRESHOLD"] > configuration["MAX_THRESHOLD"]:
        raise ValueError("MIN_THRESHOLD cannot exceed MAX_THRESHOLD")
    if configuration["MASK_BEAM_CENTER_RADIUS"] < 0.0:
        raise ValueError("MASK_BEAM_CENTER_RADIUS must be non-negative")
    return configuration


def load_configuration(path: str | Path) -> dict[str, Any]:
    """Load, normalize, and validate a TOML configuration file.

    Parameters
    ----------
    path
        Path to the TOML configuration file.

    Returns
    -------
    dict[str, Any]
        A normalized, complete configuration.

    Raises
    ------
    ValueError
        If the TOML file is malformed or contains invalid settings.
    """
    path = Path(path)
    try:
        with path.open("rb") as configuration_file:
            values = tomllib.load(configuration_file)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"Invalid TOML configuration: {path}") from error
    return merge_configuration(values)


def parse_arguments(arguments: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse the command-line configuration and output paths.

    Parameters
    ----------
    arguments
        Arguments to parse. If omitted, parse the process command line.

    Returns
    -------
    argparse.Namespace
        Parsed configuration and output paths.
    """
    parser = argparse.ArgumentParser(
        prog="prepare_sensitivities_eqsans",
        description="Prepare an EQSANS sensitivity map from a TOML configuration file.",
        epilog=f"For configuration details and usage instructions, see:\n{DOCUMENTATION_URL}",
    )
    parser.add_argument("configuration", type=Path, help="TOML configuration file.")
    parser.add_argument("output", type=Path, help="Path of the sensitivity map to create.")
    return parser.parse_args(arguments)


def main(arguments: Optional[Sequence[str]] = None) -> None:
    """Prepare an EQSANS sensitivity map using a TOML configuration file.

    Parameters
    ----------
    arguments
        Command-line arguments. If omitted, parse the process command line.

    Returns
    -------
    None
        The configured sensitivity map is written to the requested output path.
    """
    arguments = parse_arguments(arguments)
    configuration = load_configuration(arguments.configuration)
    preparer = PrepareSensitivityCorrection(component="detector1")
    preparer.set_flood_runs(cast(RunInput, configuration["FLOOD_RUNS"]))
    preparer.set_direct_beam_runs(cast(RunInput, configuration["DIRECT_BEAM_RUNS"]))
    preparer.set_masks(configuration["UNIVERSAL_MASK"], configuration["MASKED_PIXELS"])
    preparer.set_beam_center_radius(configuration["MASK_BEAM_CENTER_RADIUS"])
    preparer.set_pixel_calibration_flag(configuration["PIXEL_CALIBRATION"])
    preparer.set_solid_angle_correction_flag(configuration["SOLID_ANGLE_CORRECTION"])
    preparer.set_sample_detector_offsets(configuration["SAMPLE_OFFSET"], configuration["DETECTOR_OFFSET"])
    preparer.set_tof_cuts(configuration["LOW_TOF_CLIP"], configuration["HIGH_TOF_CLIP"])

    for option_name, setter in (
        ("DARK_CURRENT_RUNS", preparer.set_dark_current_runs),
        ("BLOCKED_BEAM_RUNS", preparer.set_blocked_beam_runs),
        ("EMPTY_CELL_RUNS", preparer.set_empty_cell_runs),
    ):
        runs = configuration[option_name]
        if runs is not None:
            setter(cast(RunInput, runs))

    if any(
        configuration[name] is not None
        for name in (
            "TRANSMISSION_FLOOD_RUNS",
            "TRANSMISSION_EMPTY_CELL_RUNS",
            "TRANSMISSION_REFERENCE_RUNS",
        )
    ):
        preparer.set_transmission_correction(
            cast(RunInput, configuration["TRANSMISSION_FLOOD_RUNS"]),
            cast(RunInput, configuration["TRANSMISSION_EMPTY_CELL_RUNS"]),
            cast(RunInput, configuration["TRANSMISSION_REFERENCE_RUNS"]),
        )

    if configuration["SCALE_COMPONENTS"] is not None:
        preparer.scale_components = configuration["SCALE_COMPONENTS"]

    preparer.execute(
        configuration["MOVING_DETECTORS"],
        configuration["MIN_THRESHOLD"],
        configuration["MAX_THRESHOLD"],
        str(arguments.output),
    )


if __name__ == "__main__":
    main()
