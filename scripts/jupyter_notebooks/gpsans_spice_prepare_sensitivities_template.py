"""
    SANS sensitivities preparation script

    # goal
    1. implement a universal mask_beam_center(flood_ws, beam_center_mask=None, beam_center_ws=None)
       for 3 types of mask
    2. add option for wing/main detector for BIOSANS:w


"""
import os
import sys
import warnings
import h5py
import numpy as np
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
warnings.simplefilter(action="ignore", category=FutureWarning)

IPTS = 828
EXPERIMENT = 280

# Set the directory for already converted SPICE files
NEXUS_DIR = os.path.join(f'/HFIR/CG2/IPTS-{IPTS}/shared', f'Exp{EXPERIMENT}')
# Check
if not os.path.exists(NEXUS_DIR):
    print(f'[ERROR] Converted NeXus-SPICE dirctory {NEXUS_DIR} does not exist')

# Input Flood Runs
FLOOD_RUNS = (38, 1), (40, 1), (42, 1)  # list of 2 tuple as (scan number, pt number)

# Pixel calibration
# Pixel calibration: False/True (default database)/user specified calibration database
# PIXEL_CALIBRATION = None
PIXEL_CALIBRATION = '/HFIR/CG2/IPTS-828/shared/pixel_calibration/runs_1_111/pixel_calibration.json'


# BIO-SANS detector
WING_DETECTOR = False  # this is main detector

# About Masks
# CG3 Main:
DIRECT_BEAM_RUNS = (37, 1), (39, 1), (41, 1)  # list of 2 tuple as (scan number, pt number)
# Beam center size
MASK_BEAM_CENTER_RADIUS = 140  # mm
BEAM_CENTER_MASKS = None

# Dark current
DARK_CURRENT_RUNS = None  # No mask, no solid angle

# Transmission run
TRANSMISSION_REFERENCE_RUNS = None  # GG3 main
# Transmission flood run
TRANSMISSION_FLOOD_RUNS = None

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
MASKED_PIXELS = '1-8,249-256'
# Mask angle: must 2 values as min and max or None
MAIN_DET_MASK_ANGLE = None
WING_DET_MASK_ANGLE = None

# Corrections
SOLID_ANGLE_CORRECTION = True   # shall be on!
TRANSMISSION_CORRECTION = False
BEAM_TRAP_SIZE_FACTOR = 2   # For BIO-SANS masking angle only.
# Flag to do dependent correction with transmission correction
THETA_DEPENDENT_CORRECTION = True

# If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 1.5

# Output
FILE_SURFIX = f'spice'

# --------------  END OF USER INPUTS --------------

# Determine sensitivities file name
if PIXEL_CALIBRATION is None:
    FILE_SURFIX += '_nobar'
else:
    FILE_SURFIX += '_bar'
SENSITIVITY_FILE = os.path.join('/HFIR/CG2/shared/drt_sensitivity/', f'sens_{INSTRUMENT}_{FILE_SURFIX}.nxs')

# --------------  DO NOT CHANGE ANY CODE BELOW THIS LINE.  THANKS! --------------------------

# Load data files
INSTRUMENT = 'CG2'
if INSTRUMENT != 'CG2':
    print('Instrument {} is not supported.  Supported are {}'
          ''.format(INSTRUMENT, 'CG2, EQSANS, CG3'))
    sys.exit(-1)

preparer = PrepareSensitivityCorrection(INSTRUMENT, WING_DETECTOR)
# Load flood runs
# map the run number to file name as it is SPICE
flood_runs = [os.path.join(NEXUS_DIR, f'CG2_{EXPERIMENT:04}{scan:04}{pt:04}.nxs.h5') for scan, pt in FLOOD_RUNS]
preparer.set_flood_runs(flood_runs)

# Process beam center runs
if DIRECT_BEAM_RUNS is not None:
    transmission_runs = [os.path.join(NEXUS_DIR, f'CG2_{EXPERIMENT:04}{scan:04}{pt:04}.nxs.h5')
                         for scan, pt in DIRECT_BEAM_RUNS]
    preparer.set_direct_beam_runs(transmission_runs)

# Set extra masks
preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS,
                   wing_det_mask_angle=WING_DET_MASK_ANGLE,
                   main_det_mask_angle=MAIN_DET_MASK_ANGLE)

# Set beam center radius
if MASK_BEAM_CENTER_RADIUS is not None:
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)
else:
    raise RuntimeError('MASK BEAM CENTER RADIUS must be set')

# Transmission
if TRANSMISSION_REFERENCE_RUNS is not None:
    raise RuntimeError('Transmission correction has not been adapted.')
    preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                         transmission_reference_run=TRANSMISSION_REFERENCE_RUNS,
                                         beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

# Dark runs
if DARK_CURRENT_RUNS is not None:
    preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

# Pixel calibration
if PIXEL_CALIBRATION:
    print(f'Pixel calibration: {PIXEL_CALIBRATION}')
    preparer.set_pixel_calibration_flag(PIXEL_CALIBRATION)

# Solid angle
preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

# Run: since it is for SPICE file, it is enforced to use IDF from NeXus
preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
                 enforce_use_nexus_idf=True)

# Information
print(f'Generated sensitivity file: {SENSITIVITY_FILE}')
# Load and print out some information
with h5py.File(SENSITIVITY_FILE) as sens:
    sens_values = sens['mantid_workspace_1']['workspace']['values'][()]
    print(f'Number of NaNs = {len(np.where(np.isnan(sens_values))[0])}')
