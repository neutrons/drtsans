import os
import h5py
import numpy as np
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection
from drtsans.mono.spice_data import SpiceRun
from typing import List, Optional, Union

MY_BEAM_LINE = "CG2"


def prepare_spice_sensitivities_correction(
    flood_spice_runs: List[SpiceRun],
    direct_beam_spice_runs: Union[List[SpiceRun], None],
    moving_detectors_methods: bool,
    min_count_threshold: float,
    max_count_threshold: float,
    nexus_dir: Union[str, None],
    mask_file: Union[str, None],
    masked_pixels: Union[str, None],
    beam_center_mask_radius: float,
    output_dir: Union[str, None],
    file_suffix: str = "spice",
    pixel_calibration_file: Union[str, None] = None,
    scale_components: Optional[dict] = None,
    solid_angle_correction: bool = True,
) -> str:
    """

    Parameters
    ----------
    flood_spice_runs: ~list
        list of SpiceRun
    direct_beam_spice_runs: ~list or None
        list of direct beam run (i.e., transmission run) in form of SpiceRun
    moving_detectors_methods: bool
        flag of sensitivity preparation method: if True, it is with moving detector. otherwise, detector patch
    min_count_threshold: float
        minimum threshold of allowed (normalized) counts
    max_count_threshold
    nexus_dir
    mask_file: str, None
        mask file applied to data
    masked_pixels: str, None
        set of pixels to mask
    beam_center_mask_radius: float
        radius in unit (mm) of beam center to mask for transmission calculation
    output_dir: str, None
        output directory, None for default as  /HFIR/{MY_BEAM_LINE}/shared/drt_sensitivity/
    file_suffix:
    pixel_calibration_file: str or None
        if it is specified as a pixel calibration, include pixel calibration in the computation
    scale_components: Optional[dict]
        Dictionary of component names and scaling factors in the form of a three-element list,
        indicating rescaling of the pixels in the component along the X, Y, and Z axes.
        For instance, ``{"detector1": [1.0, 2.0, 1.0]}`` scales pixels along the Y-axis by a factor of 2,
        leaving the other pixel dimensions unchanged.
    solid_angle_correction: bool
        do solid angle correction

    """
    # Determine output directory if default
    if output_dir is None:
        output_dir = f"/HFIR/{MY_BEAM_LINE}/shared/drt_sensitivity/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Determine sensitivities file name
    if pixel_calibration_file is None:
        file_suffix += "_nobar"
    else:
        file_suffix += "_bar"
    sens_file_name = os.path.join(output_dir, f"sens_gpsans_{file_suffix}.nxs")

    # Create the sensitivity preparation correction workflow object
    preparer = PrepareSensitivityCorrection(MY_BEAM_LINE)

    # Load flood runs
    # map the run number to file name as it is SPICE
    flood_nexus_files = [spice_run.unique_nexus_name(nexus_dir, True) for spice_run in flood_spice_runs]
    preparer.set_flood_runs(flood_nexus_files)

    # Process beam center/transmission runs
    if direct_beam_spice_runs is not None:
        transmission_nexus_files = [
            spice_run.unique_nexus_name(nexus_dir, True) for spice_run in direct_beam_spice_runs
        ]
        preparer.set_direct_beam_runs(transmission_nexus_files)

    # Set extra masks
    preparer.set_masks(mask_file, masked_pixels)

    # Set beam center radius
    if beam_center_mask_radius is not None:
        preparer.set_beam_center_radius(beam_center_mask_radius)
    else:
        raise RuntimeError("MASK BEAM CENTER RADIUS must be set")

    # Pixel calibration
    if pixel_calibration_file:
        print(f"Pixel calibration: {pixel_calibration_file}")
        preparer.set_pixel_calibration_flag(pixel_calibration_file)

    # Component scaling
    if scale_components:
        preparer.scale_components = scale_components

    # Solid angle
    preparer.set_solid_angle_correction_flag(solid_angle_correction)

    # Run: since it is for SPICE file, it is enforced to use IDF from NeXus
    try:
        preparer.execute(
            moving_detectors_methods,
            min_count_threshold,
            max_count_threshold,
            sens_file_name,
            enforce_use_nexus_idf=True,
        )
    except FileNotFoundError as file_error:
        raise file_error

    # Information
    print(f"Generated sensitivity file: {sens_file_name}")
    # Load and print out some information
    with h5py.File(sens_file_name) as sens:
        sens_values = sens["mantid_workspace_1"]["workspace"]["values"][()]
        print(f"Number of NaNs = {len(np.where(np.isnan(sens_values))[0])}")

    return sens_file_name
