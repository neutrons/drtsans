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


# INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'
INSTRUMENT = 'CG3'   # Main

# Input Flood Runs
# CG2: FLOOD_RUNS = 7116, 7118, 7120  # Single value integer or a list or tuple

# CG3: Main
FLOOD_RUNS = 4829
# BIO-SANS detector
WING_DETECTOR = False  # this is main detector

# About Masks
# CG3 Main:
DIRECT_BEAM_RUNS = 4827
# Beam center size
MASK_BEAM_CENTER_RADIUS = 65  # mm
BEAM_CENTER_MASKS = None

# Transmission run
TRANSMISSION_RUNS = 4828  # GG3 main
# Transmission flood run
TRANSMISSION_FLOOD_RUNS = 4829

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
# CG2: MASKED_PIXELS = '1-8,249-256'
# CG3:
MASKED_PIXELS = '1-18,239-256'  # CG3
# Mask angle: must 2 values as min and max or None
MAIN_DET_MASK_ANGLE = 1.5
WING_DET_MASK_ANGLE = 57.05

# Corrections
SOLID_ANGLE_CORRECTION = False
TRANSMISSION_CORRECTION = True
# Flag to do dependent correction with transmission correction
THETA_DEPENDENT_CORRECTION = True

# If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0

# Output
FILE_SURFIX = 'm7p0'
SENSITIVITY_FILE = '/HFIR/{}/shared/sens_f{}.nxs'.format(INSTRUMENT, FILE_SURFIX)

# --------------  END OF USER INPUTS --------------

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
preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS,
                   wing_det_mask_angle=WING_DET_MASK_ANGLE,
                   main_det_mask_angle=MAIN_DET_MASK_ANGLE)

# Set beam center radius
if MASK_BEAM_CENTER_RADIUS is not None:
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

# Transmission
if TRANSMISSION_RUNS is not None:
    preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                         transmission_beam_run=TRANSMISSION_RUNS)
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

# Run
preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE)

