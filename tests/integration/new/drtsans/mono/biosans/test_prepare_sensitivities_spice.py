import pytest
import os
import numpy as np
import warnings
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from drtsans.mono.spice_data import SpiceRun
warnings.simplefilter(action="ignore", category=FutureWarning)


def test_main_detector(reference_dir):
    """
        SANS sensitivities preparation script from SPICE file

        Test case for CG3 main detector

    Flood for main detector at 7m and wing detector at 12.2°-
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0009_0001.xml

    Empty Beam for Transmission Reference at 7m and 12.2° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0010_0001.xml

    Flood for wing detector at 1.4° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0020_0001.xml

    Empty Beam for Transmission Reference for wing detector at 1.4° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0016_0001.xml

    Dark Current for all configurations above -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0022_0001.xml
    """
    output_dir = os.getcwd()

    CG3 = 'CG3'  # Main

    # IPTS
    IPTS = 17241

    EXPERIMENT = 549

    # CG3: Main
    FLOOD_RUN = (9, 1)  # BIO-SANS detector
    WING_DETECTOR = False  # this is main detector

    # About Masks
    # CG3 Main beam center file/Empty beam.  It is allowed to be left blank
    DIRECT_BEAM_RUN = None
    # Beam center size
    MASK_BEAM_CENTER_RADIUS = 65  # mm
    # BEAM_CENTER_MASKS = None

    # Dark current
    DARK_CURRENT_RUN = (22, 1)

    # Transmission empty beam
    OPEN_BEAM_TRANSMISSION = (10, 1)
    # Transmission flood run
    TRANSMISSION_FLOOD_RUN = FLOOD_RUN

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

    # # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    # MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Output
    FILE_SURFIX = 'wing' if WING_DETECTOR else 'main'
    SENSITIVITY_FILE = os.path.join(output_dir, f'{CG3}_sens_{FILE_SURFIX}{FLOOD_RUN}sac_tdc7m.nxs')

    # --------------  END OF USER INPUTS --------------

    # Convert SPICE file to NeXus file
    flood_run = SpiceRun(CG3, IPTS, EXPERIMENT, FLOOD_RUN[0], FLOOD_RUN[1])
    direct_beam_run = SpiceRun(CG3, IPTS, EXPERIMENT, DIRECT_BEAM_RUN[0], DIRECT_BEAM_RUN[1]) if DIRECT_BEAM_RUN \
        else None
    open_beam_transmission = SpiceRun(CG3, IPTS, EXPERIMENT, OPEN_BEAM_TRANSMISSION[0],
                                      OPEN_BEAM_TRANSMISSION[1]) if OPEN_BEAM_TRANSMISSION else None
    transmission_flood_run = SpiceRun(CG3, IPTS, EXPERIMENT, TRANSMISSION_FLOOD_RUN[0], TRANSMISSION_FLOOD_RUN[1])
    dark_current_run = SpiceRun(CG3, IPTS, EXPERIMENT, DARK_CURRENT_RUN[0],
                                DARK_CURRENT_RUN[1]) if DARK_CURRENT_RUN else None

    prepare_spice_sensitivities_correction(WING_DETECTOR, flood_run,
                                           direct_beam_run, dark_current_run,
                                           SOLID_ANGLE_CORRECTION,
                                           transmission_flood_run, open_beam_transmission, BEAM_TRAP_SIZE_FACTOR,
                                           THETA_DEPENDENT_CORRECTION,
                                           UNIVERSAL_MASK, MASKED_PIXELS, MASK_BEAM_CENTER_RADIUS,
                                           MAIN_DET_MASK_ANGLE, WING_DET_MASK_ANGLE,
                                           MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
                                           nexus_dir=reference_dir.new.biosans)

    # Verify
    gold_sens_file = os.path.join(reference_dir.new.biosans, 'CG3_sens_main_exp549_scan9.nxs')
    assert os.path.exists(gold_sens_file)
    verify_results(SENSITIVITY_FILE, gold_sens_file)


def prepare_spice_sensitivities_correction(WING_DETECTOR, flood_run,
                                           direct_beam_run, dark_current_run,
                                           SOLID_ANGLE_CORRECTION,
                                           transmission_flood_run, open_beam_transmission, BEAM_TRAP_SIZE_FACTOR,
                                           THETA_DEPENDENT_CORRECTION,
                                           UNIVERSAL_MASK, MASKED_PIXELS, MASK_BEAM_CENTER_RADIUS,
                                           MAIN_DET_MASK_ANGLE, WING_DET_MASK_ANGLE,
                                           MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
                                           nexus_dir=None):

    CG3 = 'CG3'

    # Set up sensitivities preparation configurations
    preparer = PrepareSensitivityCorrection(CG3, WING_DETECTOR)
    # Load flood runs
    preparer.set_flood_runs([flood_run.unique_nexus_name(nexus_dir, True)])

    # Process beam center runs
    if direct_beam_run is not None:
        preparer.set_direct_beam_runs([direct_beam_run.unique_nexus_name(None, True)])

    # Set extra masks
    preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS,
                       wing_det_mask_angle=WING_DET_MASK_ANGLE,
                       main_det_mask_angle=MAIN_DET_MASK_ANGLE)

    # Set beam center radius
    if MASK_BEAM_CENTER_RADIUS is not None:
        preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)

    # Transmission
    if open_beam_transmission is not None:
        trans_flood_file = transmission_flood_run.unique_nexus_name(nexus_dir, True)
        trans_ref_file = open_beam_transmission.unique_nexus_name(nexus_dir, True)
        preparer.set_transmission_correction(transmission_flood_runs=[trans_flood_file],
                                             transmission_reference_runs=[trans_ref_file],
                                             beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
        preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)

    # Dark runs
    if dark_current_run is not None:
        dark_curr_file = dark_current_run.unique_nexus_name(nexus_dir, True)
        preparer.set_dark_current_runs([dark_curr_file])

    preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)

    # Run
    MOVING_DETECTORS = False
    preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
                     enforce_use_nexus_idf=True, debug_mode=True)

    return


def test_wing_detector(reference_dir):
    """
        SANS sensitivities preparation script from SPICE file

        Test case for CG3 main detector

    Flood for wing detector at 1.4° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0020_0001.xml

    Empty Beam for Transmission Reference for wing detector at 1.4° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0016_0001.xml

    Dark Current for all configurations above -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0022_0001.xml
    """
    output_dir = os.getcwd()

    CG3 = 'CG3'  # Main

    # IPTS
    IPTS = 17241

    EXPERIMENT = 549

    # CG3: Main
    FLOOD_RUN = (20, 1)  # BIO-SANS detector
    WING_DETECTOR = True  # this is main detector

    # About Masks
    # CG3 Main beam center file/Empty beam.  It is allowed to be left blank
    DIRECT_BEAM_RUN = None
    # Beam center size
    MASK_BEAM_CENTER_RADIUS = 65  # mm
    # BEAM_CENTER_MASKS = None


    # Transmission empty beam
    OPEN_BEAM_TRANSMISSION = (16, 1)
    # Transmission flood run
    TRANSMISSION_FLOOD_RUN = FLOOD_RUN

    # Dark current
    DARK_CURRENT_RUN = (22, 1)

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

    # # If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
    # MOVING_DETECTORS = False

    # THRESHOLD
    MIN_THRESHOLD = 0.5
    MAX_THRESHOLD = 2.0

    # Output
    FILE_SURFIX = 'wing' if WING_DETECTOR else 'main'
    SENSITIVITY_FILE = os.path.join(output_dir, f'{CG3}_sens_{FILE_SURFIX}{FLOOD_RUN}sac_tdc7m.nxs')

    # --------------  END OF USER INPUTS --------------

    # Convert SPICE file to NeXus file
    flood_run = SpiceRun(CG3, IPTS, EXPERIMENT, FLOOD_RUN[0], FLOOD_RUN[1])
    direct_beam_run = SpiceRun(CG3, IPTS, EXPERIMENT, DIRECT_BEAM_RUN[0], DIRECT_BEAM_RUN[1]) if DIRECT_BEAM_RUN \
        else None
    open_beam_transmission = SpiceRun(CG3, IPTS, EXPERIMENT, OPEN_BEAM_TRANSMISSION[0],
                                      OPEN_BEAM_TRANSMISSION[1]) if OPEN_BEAM_TRANSMISSION else None
    transmission_flood_run = SpiceRun(CG3, IPTS, EXPERIMENT, TRANSMISSION_FLOOD_RUN[0], TRANSMISSION_FLOOD_RUN[1])
    dark_current_run = SpiceRun(CG3, IPTS, EXPERIMENT, DARK_CURRENT_RUN[0],
                                DARK_CURRENT_RUN[1]) if DARK_CURRENT_RUN else None

    # # Set up sensitivities preparation configurations
    # preparer = PrepareSensitivityCorrection(CG3, WING_DETECTOR)
    # # Load flood runs
    # preparer.set_flood_runs([flood_run.unique_nexus_name(None, True)])
    #
    # # Process beam center runs
    # if direct_beam_run is not None:
    #     preparer.set_direct_beam_runs([direct_beam_run.unique_nexus_name(None, True)])
    #
    # # Set extra masks
    # preparer.set_masks(UNIVERSAL_MASK, MASKED_PIXELS,
    #                    wing_det_mask_angle=WING_DET_MASK_ANGLE,
    #                    main_det_mask_angle=MAIN_DET_MASK_ANGLE)
    #
    # # Set beam center radius
    # if MASK_BEAM_CENTER_RADIUS is not None:
    #     preparer.set_beam_center_radius(MASK_BEAM_CENTER_RADIUS)
    #
    # # Transmission
    # if open_beam_transmission is not None:
    #     preparer.set_transmission_correction(transmission_flood_runs=[transmission_flood_run.unique_nexus_name(None,
    #                                                                                                            True)],
    #                                          transmission_reference_runs=open_beam_transmission.unique_nexus_name(None,
    #                                                                                                               True),
    #                                          beam_trap_factor=BEAM_TRAP_SIZE_FACTOR)
    #     preparer.set_theta_dependent_correction_flag(THETA_DEPENDENT_CORRECTION)
    #
    # # Dark runs
    # if dark_current_run is not None:
    #     preparer.set_dark_current_runs([dark_current_run.unique_nexus_name(None, True)])
    #
    # preparer.set_solid_angle_correction_flag(SOLID_ANGLE_CORRECTION)
    #
    # # Run
    # print('Flag1')
    # MOVING_DETECTORS = False
    # preparer.execute(MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
    #                  enforce_use_nexus_idf=True, debug_mode=True)

    prepare_spice_sensitivities_correction(WING_DETECTOR, flood_run,
                                           direct_beam_run, dark_current_run,
                                           SOLID_ANGLE_CORRECTION,
                                           transmission_flood_run, open_beam_transmission, BEAM_TRAP_SIZE_FACTOR,
                                           THETA_DEPENDENT_CORRECTION,
                                           UNIVERSAL_MASK, MASKED_PIXELS, MASK_BEAM_CENTER_RADIUS,
                                           MAIN_DET_MASK_ANGLE, WING_DET_MASK_ANGLE,
                                           MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE,
                                           nexus_dir=reference_dir.new.biosans)

    # Verify
    gold_sens_file = os.path.join(reference_dir.new.biosans, 'CG3_sens_wing_exp549_scan20.nxs')
    assert os.path.exists(gold_sens_file)
    verify_results(SENSITIVITY_FILE, gold_sens_file)


def verify_results(test_sensitivities_file, gold_sens_file):
    """

    Parameters
    ----------
    test_sensitivities_file
    reference_dir: tests.conftest.rett

    Returns
    -------

    """
    from mantid.simpleapi import LoadNexusProcessed

    # Get gold file
    # gold_sens_file = os.path.join(reference_dir.new.gpsans, 'calibrations/sens_CG2_spice_bar.nxs')
    if not os.path.exists(gold_sens_file):
        raise RuntimeError(f'Expected (gold) sensitivities cannot be found at {gold_sens_file}')

    # Compare sensitivities
    gold_sens_ws = LoadNexusProcessed(Filename=gold_sens_file)
    test_sens_ws = LoadNexusProcessed(Filename=test_sensitivities_file)
    np.testing.assert_allclose(test_sens_ws.extractY(), gold_sens_ws.extractY())


if __name__ == '__main__':
    pytest.main([__file__])
