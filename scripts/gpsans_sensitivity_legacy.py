"""
    SANS sensitivities preparation script

    # goal
    1. implement a universal mask_beam_center(flood_ws, beam_center_mask=None, beam_center_ws=None)
       for 3 types of mask
    2. add option for wing/main detector for BIOSANS:w


"""
import sys
import os
import warnings
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection

warnings.simplefilter(action="ignore", category=FutureWarning)


def locate_cg2_spice_nexus(ipts_number, exp_number, runs, spice_nexus_sub):
    """Convert SPICE file information to previously converted NeXus path

    Parameters
    ----------
    runs: tuple
        tuple of 2-tuple as (scan, pt)

    Returns
    -------
    ~list
        list of converted-spice Nexus file path

    """
    # Create directory for NeXus file generated
    nexus_dir = os.path.join('/HFIR/CG2/', f'IPTS-{ipts_number}')
    nexus_dir = os.path.join(nexus_dir, os.path.join('shared', spice_nexus_sub))
    nexus_dir = os.path.join(nexus_dir, f'Exp{exp_number}')
    assert os.path.exists(nexus_dir), f'Converted NeXus directory {nexus_dir} cannot be found'

    # Create standard file name
    nexus_run_list = list()
    err_msg = ''
    for scan_num, pt_num in runs:
        nexus_name = f'CG2_{exp_number:03}{scan_num:04}{pt_num:04}.nxs.h5'
        nexus_path = os.path.join(nexus_dir, nexus_name)
        if not os.path.exists(nexus_path):
            err_msg += f'Scan {scan_num} Pt {pt_num} does not have converted Nexus {nexus_path}\n'
        else:
            nexus_run_list.append(nexus_path)

    # report error and raise exception
    if len(err_msg) > 0:
        raise RuntimeError(f'IPTS {ipts_number} Exp {exp_number} do not have all the SPICE converted\n{err_msg}')

    return nexus_run_list


INSTRUMENT = 'CG2'  # Main

IPTS = 828
EXPERIMENT = 280
SPICE_NEXUS_SUBDIR = ''

# Input Flood Runs
FLOOD_RUNS = (23, 1), (31, 1), (35, 1)  # Single tuple of a list of tuples.  Each tuple is (Scan, Pt)

# CG2
WING_DETECTOR = False  # this is main detector

# About Masks
# aka beam center runs
DIRECT_BEAM_RUNS = (24, 1), (28, 1), (20, 1)
# Beam center size
MASK_BEAM_CENTER_RADIUS = 140  # mm
BEAM_CENTER_MASKS = None

# Dark current
DARK_CURRENT_RUNS = None  # No mask, no solid angle

# Transmission run
TRANSMISSION_REFERENCE_RUNS = (23, 1), (27, ), (26, 1)  # Single tuple of a list of tuples.  Each tuple is (Scan, Pt)
# Transmission flood run
TRANSMISSION_FLOOD_RUNS = None

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
# CG2: MASKED_PIXELS = '1-8,249-256'
# CG3:
MASKED_PIXELS = '1-12,239-256'  # CG3
# Mask angle: must 2 values as min and max or None
MAIN_DET_MASK_ANGLE = 1.5
WING_DET_MASK_ANGLE = 57.05

# Adjust pixel heights and widths from bar-scan and tube-width calibrations for the following data:
# - flood runs
# - beam center runs
# - transmission runs

# Pixel calibration: False/True (default database)/user specified calibration database
PIXEL_CALIBRATION = '/HFIR/CG2/IPTS-828/shared/pixel_calibration/runs_1_111/pixel_calibration.json'

# Corrections
SOLID_ANGLE_CORRECTION = True
TRANSMISSION_CORRECTION = False
BEAM_TRAP_SIZE_FACTOR = 2  # For BIO-SANS masking angle only.
# Flag to do dependent correction with transmission correction
THETA_DEPENDENT_CORRECTION = False

# If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0

# Output
FILE_SURFIX = 'c488_nobar'
SENSITIVITY_FILE = '/HFIR/{}/shared/drt_sensitivity/sens_f{}.nxs'.format(INSTRUMENT, FILE_SURFIX)

# --------------  END OF USER INPUTS --------------

# --------------  DO NOT CHANGE ANY CODE BELOW THIS LINE.  THANKS! --------------------------

# Load data files
if INSTRUMENT not in ['CG2', 'CG3', 'EQSANS']:
    print('Instrument {} is not supported.  Supported are {}'
          ''.format(INSTRUMENT, 'CG2, EQSANS, CG3'))
    sys.exit(-1)

preparer = PrepareSensitivityCorrection(INSTRUMENT, WING_DETECTOR)
# Load flood runs
preparer.set_flood_runs(FLOOD_RUNS)

# Process beam center runs
if DIRECT_BEAM_RUNS is not None:
    beam_center_runs = locate_cg2_spice_nexus(IPTS, EXPERIMENT, DIRECT_BEAM_RUNS, SPICE_NEXUS_SUBDIR)
    preparer.set_direct_beam_runs(beam_center_runs)

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
    preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                         transmission_reference_run=TRANSMISSION_REFERENCE_RUNS,
                                         beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

# Dark runs
if DARK_CURRENT_RUNS is not None:
    preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

# Pixel calibration
preparer.set_pixel_calibration_flag(PIXEL_CALIBRATION)

# Solid angle
preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

# Run
preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE)


