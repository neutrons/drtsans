"""
    Test EASANS sensitivities preparation algorithm
"""
import pytest
import numpy as np
import os
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from mantid.simpleapi import LoadNexusProcessed


def verify_sensitivities_file(test_sens_file, gold_sens_file, atol=None):
    """
    """
    if atol is None:
        atol = 3E-5

    # Load processed NeXus files from tests and gold result
    test_sens_ws = LoadNexusProcessed(Filename=test_sens_file)
    gold_sens_ws = LoadNexusProcessed(Filename=gold_sens_file)

    # Compare number of spectra
    assert test_sens_ws.getNumberHistograms() == gold_sens_ws.getNumberHistograms()

    # Verify sensitivity value
    test_y = test_sens_ws.extractY().flatten()
    gold_y = gold_sens_ws.extractY().flatten()
    np.testing.assert_allclose(gold_y, test_y, atol=atol, equal_nan=True)

    # Verify sensitivity error
    test_e = test_sens_ws.extractE().flatten()
    gold_e = gold_sens_ws.extractE().flatten()
    np.testing.assert_allclose(gold_e, test_e, atol=atol, equal_nan=True)


def test_eqsans_prepare_sensitivities():
    """Integration test on algorithm to prepare EQSANS' sensitivities

    Returns
    -------

    """
    # INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'
    INSTRUMENT = 'EQSANS'   # Main

    # Check whether the test shall be skipped
    if not os.path.exists('/SNS/EQSANS/IPTS-24648/nexus/EQSANS_111030.nxs.h5'):
        pytest.skip('Test files cannot be accessed.')

    # Input Flood Runs
    FLOOD_RUNS = 111030

    # Beam center
    DIRECT_BEAM_RUNS = 111042

    # Beam center size
    MASK_BEAM_CENTER_RADIUS = 65  # mm

    # Dark current
    DARK_CURRENT_RUNS = 108764  # No mask, no solid angle

    MASKED_PIXELS = '1-18,239-256'

    # Corrections
    SOLID_ANGLE_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    preparer = PrepareSensitivityCorrection(INSTRUMENT, False)

    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(None, MASKED_PIXELS)

    # Set beam center radius
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

    # Dark runs
    preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = 'IntegrateTest_EQSANS_Sens.nxs'
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD,
                     output_nexus_name=output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_eq_file = '/SNS/EQSANS/shared/sans-backend/data/new/ornl' \
                   '/sans/sensitivities/EQSANS_sens_patched.nxs'

    verify_sensitivities_file(output_sens_file, gold_eq_file)


def test_cg3_main_prepare_sensitivities():
    """Integration test on algorithms to prepare sensitivities for BIOSANS's main detector

    Returns
    -------

    """
    # Check whether the test shall be skipped
    if not os.path.exists('/HFIR/CG3/IPTS-23782/nexus/CG3_4829.nxs.h5'):
        pytest.skip('Test files of CG3 cannot be accessed.')

    INSTRUMENT = 'CG3'  # Main

    # CG3: Main
    FLOOD_RUNS = 4829

    # About Masks
    # CG3 Main:
    DIRECT_BEAM_RUNS = 4827

    # Transmission run
    TRANSMISSION_RUNS = 4828  # GG3 main
    # Transmission flood run
    TRANSMISSION_FLOOD_RUNS = 4829

    # Default mask to detector
    # CG3:
    MASKED_PIXELS = '1-18,239-256'  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 1.5
    WING_DET_MASK_ANGLE = 57.0
    BEAM_TRAP_SIZE_FACTOR = 2

    # Corrections
    SOLID_ANGLE_CORRECTION = True
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    preparer = PrepareSensitivityCorrection(INSTRUMENT, is_wing_detector=False)
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    if DIRECT_BEAM_RUNS is not None:
        preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(None, MASKED_PIXELS,
                       wing_det_mask_angle=WING_DET_MASK_ANGLE,
                       main_det_mask_angle=MAIN_DET_MASK_ANGLE)

    # Transmission
    preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                         transmission_reference_run=TRANSMISSION_RUNS,
                                         beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = 'IntegrateTest_CG3_Main_Sens.nxs'
    preparer.execute(False, MIN_THRESHOLD, MAX_THRESHOLD,
                     output_nexus_name=output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_eq_file = '/SNS/EQSANS/shared/sans-backend/data/new/ornl' \
                   '/sans/sensitivities/CG3_Sens_Main.nxs'

    verify_sensitivities_file(output_sens_file, gold_eq_file)


def test_cg3_wing_prepare_sensitivities():
    """Integration test on algorithms to prepare sensitivities for BIOSANS's wing detector

    Returns
    -------

    """
    # Check whether the test shall be skipped
    if not os.path.exists('/HFIR/CG3/IPTS-23782/nexus/CG3_4835.nxs.h5'):
        pytest.skip('Test files of CG3 cannot be accessed.')

    INSTRUMENT = 'CG3'  # Main

    # CG3: Wing
    FLOOD_RUNS = 4835
    # BIO-SANS detector
    WING_DETECTOR = True  # this is main detector

    # About Masks
    # CG3 Main:
    DIRECT_BEAM_RUNS = 4830

    # Transmission run
    TRANSMISSION_RUNS = 4831  # GG3 main
    # Transmission flood run
    TRANSMISSION_FLOOD_RUNS = 4835

    # CG3:
    MASKED_PIXELS = '1-18,239-256'  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 0.75
    WING_DET_MASK_ANGLE = 57.0
    BEAM_TRAP_SIZE_FACTOR = 2

    # Corrections
    SOLID_ANGLE_CORRECTION = False
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Prepare data
    preparer = PrepareSensitivityCorrection(INSTRUMENT, WING_DETECTOR)
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(None, MASKED_PIXELS,
                       wing_det_mask_angle=WING_DET_MASK_ANGLE,
                       main_det_mask_angle=MAIN_DET_MASK_ANGLE)

    # Transmission
    preparer.set_transmission_correction(transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
                                         transmission_reference_run=TRANSMISSION_RUNS,
                                         beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = 'IntegrateTest_CG3_Wing_Sens.nxs'
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_cg2_wing_file = '/SNS/EQSANS/shared/sans-backend/data/new/ornl' \
                         '/sans/sensitivities/CG3_Sens_Wing.nxs'

    verify_sensitivities_file(output_sens_file, gold_cg2_wing_file, atol=1E-7)


def test_cg2_sensitivities():
    """Integration test on algorithms to prepare sensitivities for GPSANS's
    with moving detector method

    Returns
    -------

    """
    if not os.path.exists('/HFIR/CG2/IPTS-23801/nexus/CG2_7116.nxs.h5'):
        pytest.skip('Testing file for CG2 cannot be accessed')

    INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'

    # Input Flood Runs
    FLOOD_RUNS = 7116, 7118, 7120  # Single value integer or a list or tuple

    # About Masks
    # CG2:
    DIRECT_BEAM_RUNS = 7117, 7119, 7121
    MASK_BEAM_CENTER_RADIUS = 65  # mm

    # CG2:
    MASKED_PIXELS = '1-8,249-256'

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = True

    SOLID_ANGLE_CORRECTION = True

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    preparer = PrepareSensitivityCorrection(INSTRUMENT)
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(None, MASKED_PIXELS)

    # Set beam center radius
    preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = 'IntegrateTest_CG2_MovingDet.nxs'
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_gp_file = '/SNS/EQSANS/shared/sans-backend/data/new/ornl' \
                   '/sans/sensitivities/CG2_Sens_Moving_Dets.nxs'

    verify_sensitivities_file(output_sens_file, gold_gp_file, atol=1E-7)


if __name__ == '__main__':
    pytest.main([__file__])