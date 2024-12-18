"""
Sensitivities preparation script for GP-SANS (CG2)
"""

from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection


INSTRUMENT = "CG2"

# Input Flood Runs
FLOOD_RUNS = 90091, 90093, 90095  # Single value integer or a list or tuple

# Input Direct Beam Runs
DIRECT_BEAM_RUNS = 90092, 90094, 90096

# Input Dark Current Runs
DARK_CURRENT_RUNS = None  # No mask, no solid angle

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
MASKED_PIXELS = "1-12,239-256"
# Beam center size
MASK_BEAM_CENTER_RADIUS = 140  # mm

# Adjust pixel heights and widths from bar-scan and tube-width calibrations for the following data:
# - flood runs
# - beam center runs
# - transmission runs
PIXEL_CALIBRATION = True

# Geometry Corrections
SOLID_ANGLE_CORRECTION = True

# Single or multiple flood measurements
MOVING_DETECTORS = True

# Pixel Sensitivity thresholds
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0

# Output
FILE_SUFFIX = "c504_nobar"
SENSITIVITY_FILE = f"/HFIR/{INSTRUMENT}/shared/sens_f{FILE_SUFFIX}.nxs"

# --------------  END OF USER INPUTS --------------

# --------------  DO NOT CHANGE ANY CODE BELOW THIS LINE.  THANKS! --------------------------

# Load data files
preparer = PrepareSensitivityCorrection(instrument=INSTRUMENT, component="detector1")

# Load flood runs
preparer.set_flood_runs(FLOOD_RUNS)

# Process beam center runs
if DIRECT_BEAM_RUNS is not None:
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

# Set extra masks
preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS)

# Set beam center radius
if MASK_BEAM_CENTER_RADIUS is not None:
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)
else:
    raise RuntimeError("MASK BEAM CENTER RADIUS must be set")

# Dark runs
if DARK_CURRENT_RUNS is not None:
    preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

# Pixel calibration
preparer.set_pixel_calibration_flag(PIXEL_CALIBRATION)

# Solid angle
preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

# Run
preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE)
