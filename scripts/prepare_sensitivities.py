"""
    SANS sensitivities preparation script

    # goal
    1. implement a universal mask_beam_center(flood_ws, beam_center_mask=None, beam_center_ws=None)
       for 3 types of mask
    2. add option for wing/main detector for BIOSANS:w


"""
import sys
import warnings
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
warnings.simplefilter(action="ignore", category=FutureWarning)


INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'

# Input Flood Runs
# CG2:
FLOOD_RUNS = 7116, 7118, 7120  # Single value integer or a list or tuple
# CG3: FLOOD_RUNS = 965  # 821  # CG3

# Output
SENSITIVITY_FILE = '/HFIR/CG2/shared/sensitivity1697.nxs'

# About Masks
# CG2:
DIRECT_BEAM_RUNS = 7117, 7119, 7121
# CG3: DIRECT_BEAM_RUNS = 965  # 821
MASK_BEAM_CENTER_RADIUS = 65  # mm
BEAM_CENTER_MASKS = None

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
# CG2:
MASKED_PIXELS = '1-8,249-256'
# CG3: MASKED_PIXELS = '1-18,239-256'  # CG3
# Mask angle: must 2 values as min and max or None
MASK_ANGLES = 1.5, 100.  # 1.5, 57.0   # None

# If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0

# BIO-SANS detector
WING_DETECTOR = True

# END OF USER INPUTS

# --------------  DO NOT CHANGE ANY CODE BELOW THIS LINE.  THANKS! --------------------------

# Load data files
if INSTRUMENT not in ['CG2', 'CG3', 'EQSANS']:
    print('Instrument {} is not supported.  Supported are {}'
          ''.format(INSTRUMENT, 'CG2, EQSANS, CG3'))
    sys.exit(-1)

preparer = PrepareSensitivityCorrection(INSTRUMENT)
# Load flood runs
preparer.set_flood_runs(FLOOD_RUNS)

# Process beam center runs
if DIRECT_BEAM_RUNS is not None:
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

# Set extra masks
preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS, MASK_ANGLES)

# Set beam center radius
if MASK_BEAM_CENTER_RADIUS is not None:
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

# Run
preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE)

