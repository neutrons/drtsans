import pytest
import sys
import warnings
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
warnings.simplefilter(action="ignore", category=FutureWarning)


def test_main_detector():
    """
        SANS sensitivities preparation script from SPICE file

        Test case for CG3 main detector

    """
    CG3 = 'CG3'  # Main

    # CG3: Main
    FLOOD_RUNS = 5904  # BIO-SANS detector
    WING_DETECTOR = False  # this is main detector

    # About Masks
    # CG3 Main beam center file
    DIRECT_BEAM_RUNS = 5896
    # Beam center size
    MASK_BEAM_CENTER_RADIUS = 65  # mm
    BEAM_CENTER_MASKS = None

    # Dark current
    DARK_CURRENT_RUNS = 5884

    # Transmission empty beam
    OPEN_BEAM_TRANSMISSION = 5897
    # Transmission flood run
    TRANSMISSION_FLOOD_RUNS = FLOOD_RUNS

    # Default mask to detector
    UNIVERSAL_MASK = None  # 'Mask.XML'
    # CG3:
    MASKED_PIXELS = '1-18,249-256'  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 2.0  # 0.75# 1.5 #0.75
    WING_DET_MASK_ANGLE = 57.0
    BEAM_TRAP_SIZE_FACTOR = 1.0  # 2

    # Corrections
    SOLID_ANGLE_CORRECTION = True
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Output
    FILE_SURFIX = 'main'
    SENSITIVITY_FILE = '/SNS/users/q3n/reduce2020/sensitivity487/{}_sens_{}{}sac_tdc7m.nxs'.format(INSTRUMENT,
                                                                                                   FILE_SURFIX,
                                                                                                   FLOOD_RUNS)

    # --------------  END OF USER INPUTS --------------

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
        preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS,
                       wing_det_mask_angle=WING_DET_MASK_ANGLE,
                       main_det_mask_angle=MAIN_DET_MASK_ANGLE)

    # Set beam center radius
    if MASK_BEAM_CENTER_RADIUS is not None:
        preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

    # Transmission
    if OPEN_BEAM_TRANSMISSION is not None:
        preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                             transmission_reference_run=OPEN_BEAM_TRANSMISSION,
                                             beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
        preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    # Dark runs
    if DARK_CURRENT_RUNS is not None:
        preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE)


if __name__ == '__main__':
    pytest.main([__file__])