"""
    Test EASANS sensitivities preparation algorithm
"""
import pytest
import numpy as np
import os
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from mantid.simpleapi import LoadNexusProcessed


def verify_sensitivities_file(test_sens_file, gold_sens_file):
    """
    """
    # Load processed NeXus files from tests and gold result
    test_sens_ws = LoadNexusProcessed(Filename=test_sens_file)
    gold_sens_ws = LoadNexusProcessed(Filename=gold_sens_file)

    # Compare number of spectra
    assert test_sens_ws.getNumberHistograms() == gold_sens_ws.getNumberHistograms()

    # Verify sensitivity value
    test_y = test_sens_ws.extractY().flatten()
    gold_y = gold_sens_ws.extractY().flatten()
    np.testing.assert_allclose(gold_y, test_y, atol=3E-4, equal_nan=True)

    # Verify sensitivity error
    test_e = test_sens_ws.extractE().flatten()
    gold_e = gold_sens_ws.extractE().flatten()
    np.testing.assert_allclose(gold_e, test_e, atol=3E-5, equal_nan=True)


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

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = 'IntegrateTest_EQSANS_Sens.nxs'
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD,
                     output_nexus_name=output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_eq_file = '/SNS/snfs1/instruments/EQSANS/shared/sans-backend/data/new/ornl' \
                   '/sans/sensitivities/EQSANS_sens_patched.nxs'

    verify_sensitivities_file(output_sens_file, gold_eq_file)


if __name__ == '__main__':
    pytest.main([__file__])
