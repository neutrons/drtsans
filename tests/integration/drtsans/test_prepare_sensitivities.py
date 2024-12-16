"""
    Test EASANS sensitivities preparation algorithm
"""

import pytest
from unittest.mock import patch as mock_patch
import numpy as np
import os
from os.path import join as path_join
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from drtsans.mono.biosans.prepare_sensitivities_correction import (
    PrepareSensitivityCorrection as PrepareSensitivityCorrectionBiosans,
)
from mantid.kernel import amend_config
from mantid.simpleapi import LoadNexusProcessed
from mantid.simpleapi import DeleteWorkspace
from tempfile import mktemp


def _mock_LoadEventNexus(*args, **kwargs):
    # Substitute LoadEventNexus with LoadNexusProcessed because our synthetic files were created with SaveNexus
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


def _mock_LoadEventAsWorkspace2D(*args, **kwargs):
    # Similarly for LoadEventAsWorkspace2D
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


def verify_sensitivities_file(test_sens_file, gold_sens_file):
    """ """

    # Load processed NeXus files from tests and gold result
    test_sens_ws = LoadNexusProcessed(Filename=test_sens_file)
    gold_sens_ws = LoadNexusProcessed(Filename=gold_sens_file)

    # Compare number of spectra
    assert test_sens_ws.getNumberHistograms() == gold_sens_ws.getNumberHistograms()

    # Get sensitivity value
    test_y = test_sens_ws.extractY().flatten()
    gold_y = gold_sens_ws.extractY().flatten()

    # Get sensitivity error
    test_e = test_sens_ws.extractE().flatten()
    gold_e = gold_sens_ws.extractE().flatten()

    # Verify sensitivity value
    abs_diff = np.abs(test_y - gold_y)
    threshold = (test_e + gold_e) / 2
    ratio = abs_diff / threshold
    finite_filtered_ratio = ratio[np.isfinite(ratio)]
    np.testing.assert_array_less(finite_filtered_ratio, 0.1)

    # Verify sensitivity error
    np.testing.assert_allclose(gold_e, test_e, atol=1e-2, equal_nan=True)


@pytest.mark.mount_eqsans
def test_eqsans_prepare_sensitivities(has_sns_mount, reference_dir, cleanfile):
    """Integration test on algorithm to prepare EQSANS' sensitivities

    Returns
    -------

    """
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # INSTRUMENT = 'CG2'  # 'CG2'  # From 'EQSANS', 'CG3'
    INSTRUMENT = "EQSANS"  # Main

    # Input Flood Runs
    FLOOD_RUNS = os.path.join(reference_dir.eqsans, "EQSANS_111030.nxs.h5")

    # Beam center
    DIRECT_BEAM_RUNS = os.path.join(reference_dir.eqsans, "EQSANS_111042.nxs.h5")  # 111042

    # Beam center size
    MASK_BEAM_CENTER_RADIUS = 65  # mm

    # Dark current: No mask, no solid angle
    DARK_CURRENT_RUNS = os.path.join(reference_dir.eqsans, "EQSANS_108764.nxs.h5")  # 108764

    MASKED_PIXELS = "1-18,239-256"

    # Corrections
    SOLID_ANGLE_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

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

    # Dark runs
    preparer.set_dark_current_runs(DARK_CURRENT_RUNS)

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    # Absolute path overrides saving to the default output directory selected by the developer in Mantid's preferences.
    output_sens_file = mktemp(suffix="nxs", prefix="meta_overwrite_test1")
    print("[DEBUG] Output file: {}".format(output_sens_file))
    cleanfile(output_sens_file)
    # output_sens_file = '/tmp/IntegrateTest_EQSANS_Sens.nxs'
    preparer.execute(
        MOVING_DETECTORS,
        MIN_THRESHOLD,
        MAX_THRESHOLD,
        output_nexus_name=output_sens_file,
    )

    # Verify file existence
    assert os.path.exists(output_sens_file), "Output sensitivity file {} cannot be found".format(output_sens_file)

    # Verify value
    gold_eq_file = os.path.join(reference_dir.sans, "sensitivities", "EQSANS_sens_patched_20200602.nxs")

    verify_sensitivities_file(output_sens_file, gold_eq_file)

    # Clean
    os.remove(output_sens_file)

    # NOTE:
    # mysterious leftover workspace from this test
    # BC_EQSANS_/SNS/EQSANS/shared/sans-backend/data/ornl/sans/sns/eqsans/EQSANS_111042.nxs.h5:	37.114117 MB
    # EQSANS_111030:	43.589589 MB
    # EQSANS_111030_sensitivity:	1.179928 MB
    # EQSANS_111030_sensitivity_new:	22.355925 MB
    # gold_sens_ws:	22.355492 MB
    # test_sens_ws:	22.355925 MB
    DeleteWorkspace("BC_EQSANS_111042")
    DeleteWorkspace("EQSANS_111030")
    DeleteWorkspace("EQSANS_111030_sensitivity")
    DeleteWorkspace("EQSANS_111030_sensitivity_new")
    DeleteWorkspace("gold_sens_ws")
    DeleteWorkspace("test_sens_ws")


@pytest.mark.mount_eqsans
def test_cg3_main_prepare_sensitivities(has_sns_mount, tmp_path, cleanfile):
    """Integration test on algorithms to prepare sensitivities for BIOSANS's main detector

    Returns
    -------

    """
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # Check whether the test shall be skipped
    if not os.path.exists("/HFIR/CG3/IPTS-23782/nexus/CG3_4829.nxs.h5"):
        pytest.skip("Test files of CG3 cannot be accessed.")

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
    MASKED_PIXELS = "1-18,239-256"  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 1.5
    BEAM_TRAP_SIZE_FACTOR = 2

    # Corrections
    SOLID_ANGLE_CORRECTION = True
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    preparer = PrepareSensitivityCorrectionBiosans(component="detector1")
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    if DIRECT_BEAM_RUNS is not None:
        preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(
        None,
        MASKED_PIXELS,
        main_det_mask_angle=MAIN_DET_MASK_ANGLE,
    )

    # Transmission
    preparer.set_transmission_correction(
        transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
        transmission_reference_runs=TRANSMISSION_RUNS,
        beam_trap_factor=BEAM_TRAP_SIZE_FACTOR,
    )
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    # Solid angle
    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = path_join(tmp_path, "IntegrateTest_CG3_Main_Sens.nxs")
    preparer.execute(False, MIN_THRESHOLD, MAX_THRESHOLD, output_nexus_name=output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_eq_file = "/SNS/EQSANS/shared/sans-backend/data/ornl/sans/sensitivities/CG3_Sens_Main.nxs"

    verify_sensitivities_file(output_sens_file, gold_eq_file)

    # Clean
    cleanfile(output_sens_file)


@pytest.mark.mount_eqsans
def test_cg3_wing_prepare_sensitivities(has_sns_mount, tmp_path):
    """Integration test on algorithms to prepare sensitivities for BIOSANS's wing detector

    Returns
    -------

    """
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # Check whether the test shall be skipped
    if not os.path.exists("/HFIR/CG3/IPTS-23782/nexus/CG3_4835.nxs.h5"):
        pytest.skip("Test files of CG3 cannot be accessed.")

    # CG3: Wing
    FLOOD_RUNS = 4835
    # BIO-SANS detector

    # About Masks
    # CG3 Main:
    DIRECT_BEAM_RUNS = 4830

    # Transmission run
    TRANSMISSION_RUNS = 4831  # GG3 main
    # Transmission flood run
    TRANSMISSION_FLOOD_RUNS = 4835

    # CG3:
    MASKED_PIXELS = "1-18,239-256"  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 0.75
    BEAM_TRAP_SIZE_FACTOR = 2

    # Corrections
    SOLID_ANGLE_CORRECTION = True
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Prepare data
    preparer = PrepareSensitivityCorrectionBiosans(component="wing_detector")
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(
        None,
        MASKED_PIXELS,
        main_det_mask_angle=MAIN_DET_MASK_ANGLE,
    )

    # Transmission
    preparer.set_transmission_correction(
        transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
        transmission_reference_runs=TRANSMISSION_RUNS,
        beam_trap_factor=BEAM_TRAP_SIZE_FACTOR,
    )
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = path_join(tmp_path, "IntegrateTest_CG3_Wing_Sens.nxs")
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_cg2_wing_file = "/SNS/EQSANS/shared/sans-backend/data/ornl/sans/sensitivities/CG3_Sens_Wing.nxs"

    verify_sensitivities_file(output_sens_file, gold_cg2_wing_file)

    # Clean
    os.remove(output_sens_file)

    # NOTE:
    # mysterious leftover workspaces in memory
    # BC_CG3_CG3_4830:	2.763785 MB
    # BIOSANS_4835: 44.614937 MB
    # BIOSANS_4835_sensitivity:	2.162968 MB
    # BIOSANS_4835_sensitivity_new:	5.686553 MB
    # gold_sens_ws:	5.686328 MB
    # test_sens_ws:	5.686553 MB
    # TRANS_CG3_4831:	2.762857 MB
    # TRANS_CG3_4835:	5.687177 MB
    DeleteWorkspace("BC_CG3_4830")
    DeleteWorkspace("BIOSANS_4835")
    DeleteWorkspace("BIOSANS_4835_sensitivity")
    DeleteWorkspace("BIOSANS_4835_sensitivity_new")
    DeleteWorkspace("gold_sens_ws")
    DeleteWorkspace("test_sens_ws")
    DeleteWorkspace("TRANS_CG3_4831")
    DeleteWorkspace("TRANS_CG3_4835")


@pytest.mark.datarepo
@mock_patch("drtsans.load.LoadEventAsWorkspace2D", new=_mock_LoadEventAsWorkspace2D)
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
def test_cg3_midrange_prepare_sensitivities(biosans_synthetic_sensitivity_dataset, tmp_path, cleanfile):
    """Integration test on algorithms to prepare sensitivities for BIOSANS's midrange detector"""
    # CG3: Mid
    FLOOD_RUNS = 4835
    # BIO-SANS detector

    # About Masks
    # CG3 Main:
    DIRECT_BEAM_RUNS = 4830

    # Transmission run
    TRANSMISSION_RUNS = 4831
    # Transmission flood run
    TRANSMISSION_FLOOD_RUNS = 4835

    # CG3:
    MASKED_PIXELS = "1-18,239-256"  # CG3
    # Mask angle: must 2 values as min and max or None
    MAIN_DET_MASK_ANGLE = 0.75
    BEAM_TRAP_SIZE_FACTOR = 2

    # Corrections
    SOLID_ANGLE_CORRECTION = True
    # Flag to do dependent correction with transmission correction
    THETA_DEPENDENT_CORRECTION = True

    # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Prepare data
    preparer = PrepareSensitivityCorrectionBiosans(component="midrange_detector")
    # Load flood runs
    preparer.set_flood_runs(FLOOD_RUNS)

    # Process beam center runs
    preparer.set_direct_beam_runs(DIRECT_BEAM_RUNS)

    # Set extra masks
    preparer.set_masks(
        None,
        MASKED_PIXELS,
        main_det_mask_angle=MAIN_DET_MASK_ANGLE,
    )

    # Transmission
    preparer.set_transmission_correction(
        transmission_flood_runs=TRANSMISSION_FLOOD_RUNS,
        transmission_reference_runs=TRANSMISSION_RUNS,
        beam_trap_factor=BEAM_TRAP_SIZE_FACTOR,
    )
    preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    output_sens_file = path_join(tmp_path, "IntegrateTest_CG3_Mid_Sens.nxs")
    with amend_config(data_dir=biosans_synthetic_sensitivity_dataset["data_dir"]):
        preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)
    cleanfile(output_sens_file)

    # Verify value
    # NOTE: since we are using synthetic data, there is no gold/reference file we can use to verify
    #       that the content is correct. Once we have the correct data file (instead of mock one),
    #       we should generate the gold file and enable the checking.
    # gold_cg2_wing_file = "/SNS/EQSANS/shared/sans-backend/data/ornl" "/sans/sensitivities/CG3_Sens_Wing.nxs"

    # verify_sensitivities_file(output_sens_file, gold_cg2_wing_file, atol=1e-7)


@pytest.mark.mount_eqsans
def test_cg2_sensitivities(has_sns_mount, tmp_path):
    """Integration test on algorithms to prepare sensitivities for GPSANS's
    with moving detector method

    Returns
    -------

    """
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    if not os.path.exists("/HFIR/CG2/IPTS-23801/nexus/CG2_7116.nxs.h5"):
        pytest.skip("Testing file for CG2 cannot be accessed")

    INSTRUMENT = "CG2"  # 'CG2'  # From 'EQSANS', 'CG3'

    # Input Flood Runs
    FLOOD_RUNS = 7116, 7118, 7120  # Single value integer or a list or tuple

    # About Masks
    # CG2:
    DIRECT_BEAM_RUNS = 7117, 7119, 7121
    MASK_BEAM_CENTER_RADIUS = 65  # mm

    # CG2:
    MASKED_PIXELS = "1-8,249-256"

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
    output_sens_file = path_join(tmp_path, "IntegrateTest_CG2_MovingDet.nxs")
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, output_sens_file)

    # Verify file existence
    assert os.path.exists(output_sens_file)

    # Verify value
    gold_gp_file = "/SNS/EQSANS/shared/sans-backend/data/ornl/sans/sensitivities/CG2_Sens_Moving_Dets.nxs"

    verify_sensitivities_file(output_sens_file, gold_gp_file)

    # Clean
    os.remove(output_sens_file)

    # NOTE:
    # mysterious leftover workspaces in memory
    # BC_CG2_CG2_7117:	1.434937 MB
    # BC_CG2_CG2_7119:	1.448089 MB
    # BC_CG2_CG2_7121:	1.442393 MB
    # gold_sens_ws:	12.333224 MB
    # GPSANS_7116:	33.567113 MB
    # GPSANS_7116_processed_histo:	12.334073 MB
    # GPSANS_7118:	12.118841 MB
    # GPSANS_7118_processed_histo:	12.118841 MB
    # GPSANS_7120:	12.078681 MB
    # GPSANS_7120_processed_histo:	12.078681 MB
    # sensitivities:	1.179928 MB
    # sensitivities_new:	12.333449 MB
    # test_sens_ws:	12.333449 MB
    DeleteWorkspace("BC_CG2_7117")
    DeleteWorkspace("BC_CG2_7119")
    DeleteWorkspace("BC_CG2_7121")
    DeleteWorkspace("gold_sens_ws")
    DeleteWorkspace("GPSANS_7116")
    DeleteWorkspace("GPSANS_7116_processed_histo")
    DeleteWorkspace("GPSANS_7118")
    DeleteWorkspace("GPSANS_7118_processed_histo")
    DeleteWorkspace("GPSANS_7120")
    DeleteWorkspace("GPSANS_7120_processed_histo")
    DeleteWorkspace("sensitivities")
    DeleteWorkspace("sensitivities_new")
    DeleteWorkspace("test_sens_ws")


if __name__ == "__main__":
    pytest.main([__file__])
