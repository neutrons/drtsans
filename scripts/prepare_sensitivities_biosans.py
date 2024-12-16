"""
    Sensitivities preparation script for Bio-SANS (CG3)
"""

from drtsans.mono.biosans.prepare_sensitivities_correction import PrepareSensitivityCorrection

import os


INSTRUMENT = "CG3"
RUNCYCLE = 504
IPTS = 32383

# Input Runs, one each for the main, midrange and wing detectors
FSUF = ["m6p3", "mr2p7", "w7p25"]  # for 'FILE_SUFFIX' - indicates detector and detector position (distance/angle)
WDET = ["detector1", "midrange_detector", "wing_detector"]  # for detector
DBEAM_RUNS = [22608, 22608, 22731]  # for 'DIRECT_BEAM_RUNS' -- Beam Center Measurements
TRANS_REF_RUNS = [22609, 22609, 22732]  # for 'TRANSMISSION_REFERENCE_RUNS' -- Direct Beam Transmission Measurements
FD_RUNS = [22719, 22719, 22733]  # for 'FLOOD_RUNS' -- Water 1mm Measurements
TRANSFD_RUNS = [22719, 22719, 22733]  # for 'TRANSMISSION_FLOOD_RUNS' -- Flood Transmission Measurements
DK_RUNS = [22734, 22734, 22734]  # for 'DARK_CURRENT_RUNS' -- Dark Current
SAMPLEHOLDER_SUFFIX = (
    "URbBj"  # 'Pl' is Peltier; 'Tb' is Tumbler; 'Ti' is Titanium; 'Bj' is banjo; 'Pr' Enhanced Pressure Cell, EAP
)

# Default mask
UNIVERSAL_MASK = f"/HFIR/{INSTRUMENT}/shared/Cycle504/mask_22999.xml"  # 'Mask.XML'
MASKED_PIXELS = "1-18,239-256"  # CG3
MAIN_DET_MASK_ANGLE = 1.25  # Maximum angle for circular mask on main detector to mask direct beam spot OR NONE
BEAM_TRAP_SIZE_FACTOR = 1.0  # Defines Direct Beam pixels for transmission calculations.
BEAM_CENTER_RADIUS = 38  # mm

# Adjust pixel heights and widths from bar-scan and tube-width calibrations for the following data:
# - flood runs
# - beam center runs
# - transmission runs
PIXEL_CALIBRATION = True

# Adjust pixel positions, heights and widths Mantid algorithm ScaleInstrumentComponent
# A valid input is a dictionary of scaling triads. For instance,
#
# SCALE_COMPONENT={"detector1": [1.0, 2.0, 1.0], "wing_detector":[0.5, 1.0, 0.5]}
#
# doubles the height of pixels in "detector1" (Y-axis is the vertical axis in the detector coordinate system),
# and halves the width of pixels in "wing_detector". No changes for "midrange detector".
SCALE_COMPONENTS = None  # no change to pixel positions, heights and widths

# Geometry Corrections
SOLID_ANGLE_CORRECTION = True
# Flag to do dependent correction with transmission correction
THETA_DEPENDENT_CORRECTION = True

# Single or multiple Flood measurements
MOVING_DETECTORS = False

# Pixel Sensitivity thresholds
if PIXEL_CALIBRATION:
    MIN_THRESHOLD = 0.5  # All Pixels with Normalized counts < 0.5 are disregarded
    MAX_THRESHOLD = 1.5  # All Pixels with Normalized counts > 2.0 are disregarded
else:
    MIN_THRESHOLD = 0.5  # All Pixels with Normalized counts < 0.5 are disregarded
    MAX_THRESHOLD = 2.0  # All Pixels with Normalized counts > 2.0 are disregarded

SENSITIVITY_FILEPATH = f"/HFIR/{INSTRUMENT}/IPTS-{IPTS}/shared/RC{RUNCYCLE}/Calibration/Sens/"

# --------------  END OF USER INPUTS --------------

# --------------  DO NOT CHANGE ANY CODE BELOW THIS LINE.  THANKS! --------------------------

# Loop through MAIN, MIDRANGE & WING detector sensitivity file preparations
for (
    FILE_SUFFIX,
    DETECTOR,
    FLOOD_RUNS,
    DIRECT_BEAM_RUNS,
    TRANSMISSION_REFERENCE_RUNS,
    TRANSMISSION_FLOOD_RUNS,
    DARK_CURRENT_RUNS,
) in zip(FSUF, WDET, FD_RUNS, DBEAM_RUNS, TRANS_REF_RUNS, TRANSFD_RUNS, DK_RUNS):
    preparer = PrepareSensitivityCorrection(DETECTOR)

    # Select flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Choose beam center runs
    if DIRECT_BEAM_RUNS is not None:
        preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Reduce beam spot mask to avoid masking midrange if sensitivity of midrange is being calculated
    if DETECTOR == "midrange_detector":
        MAIN_DET_MASK_ANGLE = 0.5  # Maximum angle for circular mask on main detector to mask direct beam spot OR NONE

    # Set extra masks
    preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS, MAIN_DET_MASK_ANGLE)

    # Set beam center radius
    if preparer.beam_center_radius is None:
        preparer.set_beam_center_radius(BEAM_CENTER_RADIUS)

    # Transmission
    if TRANSMISSION_REFERENCE_RUNS is not None:
        preparer.set_transmission_correction(
            transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
            transmission_reference_runs=TRANSMISSION_REFERENCE_RUNS,
            beam_trap_factor=BEAM_TRAP_SIZE_FACTOR,
        )
        preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    # Dark runs
    if DARK_CURRENT_RUNS is not None:
        preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

    # Component scaling
    preparer.scale_components = SCALE_COMPONENTS

    # Pixel calibration
    preparer.set_pixel_calibration_flag(PIXEL_CALIBRATION)

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Output
    OUTPUT_SENSFILE = os.path.join(SENSITIVITY_FILEPATH, f"Sens_f{FLOOD_RUNS}{FILE_SUFFIX}{SAMPLEHOLDER_SUFFIX}.nxs")

    # Run
    preparer.execute(
        use_moving_detector_method=MOVING_DETECTORS,
        min_threshold=MIN_THRESHOLD,
        max_threshold=MAX_THRESHOLD,
        output_nexus_name=OUTPUT_SENSFILE,
        enforce_use_nexus_idf=True,
    )
