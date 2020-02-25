"""
    Test EASANS sensitivities preparation algorithm
"""
import pytest
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection


def test_eqsans_prepare_sensitivities():
    """Integration test on algorithm to prepare EQSANS' sensitivities

    Returns
    -------

    """
    # INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'
    INSTRUMENT = 'EQSANS'   # Main

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

    # Run
    output_sens_file = 'IntegrateTest_EQSANS_Sens.nxs'
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD,
                     output_nexus_name=output_sens_file)

if __name__ == '__main__':
    pytest.main([__file__])
