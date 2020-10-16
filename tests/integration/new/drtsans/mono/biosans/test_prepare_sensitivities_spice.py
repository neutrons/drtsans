import pytest
import os
import numpy as np
import warnings
from tempfile import mkdtemp
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from drtsans.mono.spice_data import SpiceRun
from mantid.simpleapi import LoadNexusProcessed
warnings.simplefilter(action="ignore", category=FutureWarning)


def test_main_detector(reference_dir, cleanfile):
    """Test case for CG3 main detector

    Flood for main detector at 7m and wing detector at 12.2°-
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0009_0001.xml

    Empty Beam for Transmission Reference at 7m and 12.2° -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0010_0001.xml

    Dark Current for all configurations above -
    /HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0022_0001.xml
    """
    output_dir = mkdtemp()

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


def prepare_spice_sensitivities_correction(is_wing_detector: bool,
                                           flood_run: SpiceRun,
                                           direct_beam_run: SpiceRun,
                                           dark_current_run: SpiceRun,
                                           apply_solid_angle_correction: bool,
                                           transmission_flood_run: SpiceRun,
                                           transmission_reference_run: SpiceRun,
                                           beam_trap_size_factor: float,
                                           apply_theta_dependent_correction: bool,
                                           universal_mask_file: str,
                                           pixels_to_mask: str,
                                           beam_center_mask_radius: float,
                                           main_detector_mask_angle: float,
                                           wing_detector_mask_angle: float,
                                           min_count_threshold: float,
                                           max_count_threshold: float,
                                           sensitivity_file_name: str,
                                           nexus_dir: str = None):
    """Prepare sensitivities from SPICE files
    
    Parameters
    ----------
    is_wing_detector: bool
        Flag to indicate the operation is on wing detector
    flood_run: SpiceRun
        flood run
    direct_beam_run: SpiceRun or None
        direct beam run
    dark_current_run: SpiceRun
        dark current run
    apply_solid_angle_correction: bool
        Flag to apply solid angle correction to flood run
    transmission_flood_run: SpiceRun
        transmission flood run
    transmission_reference_run: SpiceRun
        transmission reference run
    beam_trap_size_factor: float
        size factor of beam trap given by user
    apply_theta_dependent_correction: bool
        Flag to apply theta dependent correction to transmission run
    universal_mask_file: str
        path to mask file applied to all the runs
    pixels_to_mask: str
        lists of pixels (IDs) to mask
    beam_center_mask_radius: float
        radius of mask for beam center in mm
    main_detector_mask_angle: float
        angle for main detector mask
    wing_detector_mask_angle: float
        angle for wing detector mask
    min_count_threshold: float
        minimum normalized count threshold as a good pixel
    max_count_threshold: float
        maximum normalized count threshold as a good pixel
    sensitivity_file_name: str
        output file name with full path
    nexus_dir: str or None
        directory for nexus file.  None for default.

    """

    CG3 = 'CG3'

    # Set up sensitivities preparation configurations
    preparer = PrepareSensitivityCorrection(CG3, is_wing_detector)
    # Load flood runs
    preparer.set_flood_runs([flood_run.unique_nexus_name(nexus_dir, True)])

    # Process beam center runs
    if direct_beam_run is not None:
        preparer.set_direct_beam_runs([direct_beam_run.unique_nexus_name(None, True)])

    # Set extra masks
    preparer.set_masks(universal_mask_file, pixels_to_mask,
                       wing_det_mask_angle=wing_detector_mask_angle,
                       main_det_mask_angle=main_detector_mask_angle)

    # Set beam center radius
    if beam_center_mask_radius is not None:
        preparer.set_beam_center_radius(beam_center_mask_radius)

    # Transmission
    if transmission_reference_run is not None:
        trans_flood_file = transmission_flood_run.unique_nexus_name(nexus_dir, True)
        trans_ref_file = transmission_reference_run.unique_nexus_name(nexus_dir, True)
        preparer.set_transmission_correction(transmission_flood_runs=[trans_flood_file],
                                             transmission_reference_runs=[trans_ref_file],
                                             beam_trap_factor=beam_trap_size_factor)
        preparer.set_theta_dependent_correction_flag(apply_theta_dependent_correction)

    # Dark runs
    if dark_current_run is not None:
        dark_curr_file = dark_current_run.unique_nexus_name(nexus_dir, True)
        preparer.set_dark_current_runs([dark_curr_file])

    preparer.set_solid_angle_correction_flag(apply_solid_angle_correction)

    # Run
    moving_detector = False
    preparer.execute(moving_detector, min_count_threshold, max_count_threshold, sensitivity_file_name,
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


def verify_results(test_sensitivities_file: str, gold_sens_file: str):
    """Verify sensitivities of tested result from gold file
    """
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
