"""Prepare an EQSANS sensitivity correction map.

Edit the input constants for the experiment before running this script. The
script intentionally uses the public EQSANS sensitivity-preparation API so it
can serve as a starting point for production workflows.
"""

from typing import Optional, Sequence, Union, cast

from drtsans.tof.eqsans.prepare_sensitivities_correction import PrepareSensitivityCorrection


Run = Union[int, str]
RunInput = Union[Run, Sequence[Run]]
OptionalRunInput = Optional[RunInput]


# Run inputs may be one integer run number, one string run identifier/path, or a
# sequence of either. Examples: 12345, "EQSANS_12345", "/path/to/run.nxs",
# or (12345, 12346). A sequence pairs its entries with the flood runs below.
FLOOD_RUNS: RunInput = (0,)
DIRECT_BEAM_RUNS: RunInput = (0,)

# Optional auxiliary run inputs use the same forms as the flood runs. Leave
# these as None when the correction is not needed. One run is broadcast to all
# flood runs; otherwise provide one auxiliary run per flood run.
DARK_CURRENT_RUNS: OptionalRunInput = None
BLOCKED_BEAM_RUNS: OptionalRunInput = None
EMPTY_CELL_RUNS: OptionalRunInput = None

# Optional wavelength-dependent transmission runs use the same run-input forms.
# Configure all three together: transmission flood, transmission empty-cell,
# and transmission reference. Each may be one run or one run per flood run.
TRANSMISSION_FLOOD_RUNS: OptionalRunInput = None
TRANSMISSION_EMPTY_CELL_RUNS: OptionalRunInput = None
TRANSMISSION_REFERENCE_RUNS: OptionalRunInput = None

# EQSANS geometry and time-of-flight options.
SAMPLE_OFFSET = 0.0  # mm
DETECTOR_OFFSET = 0.0  # mm
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

# Sensitivity calculation and output.
# Select the sensitivity-preparation method.
# Input type: bool. Use True for moving-detector data, or False for the
# detector-patch method. Example: False.
MOVING_DETECTORS = False
# Minimum accepted normalized detector response.
# Input type: float. Pixels below this value are masked. Example: 0.5.
MIN_THRESHOLD = 0.5
# Maximum accepted normalized detector response.
# Input type: float. Pixels above this value are masked. Example: 2.0.
MAX_THRESHOLD = 2.0
OUTPUT_SENSITIVITY_FILE = "/SNS/EQSANS/shared/sensitivity/eqsans_sensitivity.nxs"


def main() -> None:
    """Prepare an EQSANS sensitivity map using the configured run inputs.

    Returns
    -------
    None
        The configured sensitivity map is written to
        ``OUTPUT_SENSITIVITY_FILE``.
    """
    preparer = PrepareSensitivityCorrection(component="detector1")
    preparer.set_flood_runs(FLOOD_RUNS)
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)
    preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS)
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)
    preparer.set_pixel_calibration_flag(PIXEL_CALIBRATION)
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)
    preparer.set_sample_detector_offsets(SAMPLE_OFFSET, DETECTOR_OFFSET)
    preparer.set_tof_cuts(LOW_TOF_CLIP, HIGH_TOF_CLIP)

    if DARK_CURRENT_RUNS is not None:
        preparer.set_dark_current_runs(DARK_CURRENT_RUNS)
    if BLOCKED_BEAM_RUNS is not None:
        preparer.set_blocked_beam_runs(BLOCKED_BEAM_RUNS)
    if EMPTY_CELL_RUNS is not None:
        preparer.set_empty_cell_runs(EMPTY_CELL_RUNS)

    if any(
        runs is not None
        for runs in (
            TRANSMISSION_FLOOD_RUNS,
            TRANSMISSION_EMPTY_CELL_RUNS,
            TRANSMISSION_REFERENCE_RUNS,
        )
    ):
        preparer.set_transmission_correction(
            cast(RunInput, TRANSMISSION_FLOOD_RUNS),
            cast(RunInput, TRANSMISSION_EMPTY_CELL_RUNS),
            cast(RunInput, TRANSMISSION_REFERENCE_RUNS),
        )

    if SCALE_COMPONENTS is not None:
        preparer.scale_components = SCALE_COMPONENTS

    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, OUTPUT_SENSITIVITY_FILE)


if __name__ == "__main__":
    main()
