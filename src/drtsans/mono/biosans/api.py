"""BIOSANS API"""

# local imports
import drtsans
from drtsans import getWedgeSelection, subtract_background, NoDataProcessedError
from drtsans.dataobjects import save_i1d, IQazimuthal
from drtsans.instruments import extract_run_number
from drtsans.iq import bin_all
from drtsans.load import move_instrument
from drtsans.mask_utils import apply_mask, load_mask
from drtsans.mono import biosans
from drtsans.mono import meta_data
from drtsans.mono.biosans import solid_angle_correction
from drtsans.mono.biosans.geometry import has_midrange_detector
from drtsans.mono.dark_current import subtract_dark_current
from drtsans.mono.load import load_events, transform_to_wavelength, set_init_uncertainties
from drtsans.mono.meta_data import get_sample_detector_offset, parse_json_meta_data, set_meta_data
from drtsans.mono.normalization import (
    normalize_by_monitor,
    normalize_by_time,
    NoMonitorMetadataError,
    ZeroMonitorCountsError,
)
from drtsans.mono.transmission import apply_transmission_correction, calculate_transmission
from drtsans.path import abspath, abspaths, allow_overwrite, registered_workspace
from drtsans.plots import plot_detector
from drtsans.samplelogs import SampleLogs
from drtsans.save_ascii import save_ascii_binned_2D
from drtsans.save_cansas import save_cansas_nx
from drtsans.sensitivity import apply_sensitivity_correction, load_sensitivity_workspace
from drtsans.settings import namedtuplefy
from drtsans.stitch import stitch_binned_profiles
from drtsans.reductionlog import savereductionlog
from drtsans.thickness_normalization import normalize_by_thickness
from drtsans.transmission import TransmissionErrorToleranceError, TransmissionNanError

# third party imports
from mantid.dataobjects import Workspace2D
from mantid.kernel import Logger
from mantid.simpleapi import (
    mtd,
    MaskDetectors,
    LoadEventNexus,
    LoadNexusProcessed,
    DeleteWorkspace,
    RemoveWorkspaceHistory,
)
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# standard imports
from collections import namedtuple
import copy
from datetime import datetime
import os
from typing import Union, List

# Functions exposed to the general user (public) API
__all__ = [
    "prepare_data",
    "load_all_files",
    "plot_reduction_output",
    "prepare_data_workspaces",
    "process_single_configuration",
    "reduce_single_configuration",
]

# silicon window to nominal position (origin) distance in meter
SI_WINDOW_NOMINAL_DISTANCE_METER = 0.071
SAMPLE_SI_META_NAME = "CG3:CS:SampleToSi"

# default transmission error tolerance
DEFAULT_TRANSMISSION_ERROR_TOLERANCE = 0.01

# setup logger
# NOTE: If logging information is not showing up, please check the mantid log level.
#       If problem persists, please visit:
#       https://docs.mantidproject.org/nightly/concepts/PropertiesFile.html#logging-properties
logger = Logger("BioSANS")


@namedtuplefy
def load_all_files(
    reduction_input: dict,
    prefix: str = "",
    load_params: dict = {},
    path: str = None,
    use_nexus_idf: bool = False,
    debug_output: bool = False,
) -> dict:
    """Load all required files at the beginning, and transform them into histograms.

    Parameters
    ----------
    reduction_input: dict
        Dictionary containing the reduction input
    prefix: str
        Prefix to be used for the workspaces
    load_params: dict
        Dictionary containing the parameters to be passed to the Load event algorithm
    path: str
        Additional path to search the for NeXus files
    use_nexus_idf: bool
        Flag to enforce to use IDF from NeXus file.  It must be true for SPICE-converted NeXus
    debug_output: bool
        Flag to save internal data for debugging

    Returns
    -------
    dict
        Dictionary containing the configuration used to load the workspaces
    """
    # a handy shortcut to the configuration parameters dictionary
    reduction_config = reduction_input["configuration"]

    instrument_name = reduction_input["instrumentName"]
    ipts = reduction_input["iptsNumber"]
    sample = reduction_input["sample"]["runNumber"]

    # on the fly check to see if mid-range detector is present in data
    reduction_input["has_midrange_detector"] = file_has_midrange_detector(
        sample=sample,
        ipts=ipts,
        instrument_name=instrument_name,
        directory=path,
    )

    # check that the configuration is consistent with the presence/absence of the midrange detector in the run
    check_overlap_stitch_configuration(reduction_input)

    sample_trans = reduction_input["sample"]["transmission"]["runNumber"]
    bkgd = reduction_input["background"]["runNumber"]
    bkgd_trans = reduction_input["background"]["transmission"]["runNumber"]
    empty = reduction_input["emptyTransmission"]["runNumber"]
    center = reduction_input["beamCenter"]["runNumber"]
    blocked_beam = reduction_config["blockedBeamRunNumber"]

    # Remove existing workspaces, this is to guarantee that all the data is loaded correctly
    # In the future this should be made optional
    ws_to_remove = [
        f"{prefix}_{instrument_name}_{run_number}_raw_histo"
        for run_number in (
            sample,
            center,
            bkgd,
            empty,
            sample_trans,
            bkgd_trans,
            blocked_beam,
        )
    ]
    ws_to_remove.append(f"{prefix}_{instrument_name}_{sample}_raw_histo_slice_group")
    ws_to_remove.append(f"{prefix}_main_sensitivity")
    ws_to_remove.append(f"{prefix}_wing_sensitivity")
    ws_to_remove.append(f"{prefix}_midrange_sensitivity")
    ws_to_remove.append(f"{prefix}_mask")
    if reduction_config["darkMainFileName"]:
        run_number = extract_run_number(reduction_config["darkMainFileName"])
        ws_to_remove.append(f"{prefix}_{instrument_name}_{run_number}_raw_histo")
    if reduction_config["darkWingFileName"]:
        run_number = extract_run_number(reduction_config["darkWingFileName"])
        ws_to_remove.append(f"{prefix}_{instrument_name}_{run_number}_raw_histo")
    if reduction_config.get("darkMidrangeFileName", None):
        run_number = extract_run_number(reduction_config["darkMidrangeFileName"])
        ws_to_remove.append(f"{prefix}_{instrument_name}_{run_number}_raw_histo")
    for ws_name in ws_to_remove:
        if registered_workspace(ws_name):
            mtd.remove(ws_name)

    # load nexus idf
    if use_nexus_idf:
        load_params["LoadNexusInstrumentXML"] = use_nexus_idf

    # Adjust pixel positions, heights and widths
    load_params["scale_components"] = reduction_config.get("scaleComponents", None)
    load_params["pixel_calibration"] = reduction_config.get("usePixelCalibration", False)

    # wave length and wave length spread
    (
        wave_length_dict,
        wave_length_spread_dict,
    ) = meta_data.parse_json_wave_length_and_spread(reduction_input)

    if reduction_config["useDefaultMask"]:
        # reduction_config["defaultMask"] is a list of python dictionaries
        default_mask = reduction_config["defaultMask"] if reduction_config["defaultMask"] is not None else []
    else:
        default_mask = []

    # check for time/log slicing
    timeslice = reduction_config["useTimeSlice"]
    logslice = reduction_config["useLogSlice"]
    if timeslice or logslice:
        if len(sample.split(",")) > 1:
            logger.error("Can't do slicing on summed data sets")
            raise ValueError("Can't do slicing on summed data sets")

    # thickness is written to sample log if it is defined...
    thickness = reduction_input["sample"]["thickness"]
    sample_aperture_diameter = reduction_config["sampleApertureSize"]  # in mm
    source_aperture_diameter = reduction_config["sourceApertureDiameter"]  # in mm

    # Parse smearing pixel size x and y for all runs
    smearing_pixel_size_x_dict = parse_json_meta_data(
        reduction_input,
        "smearingPixelSizeX",
        1e-3,
        beam_center_run=True,
        background_run=True,
        empty_transmission_run=True,
        transmission_run=True,
        background_transmission=True,
        block_beam_run=True,
        dark_current_run=True,
    )

    smearing_pixel_size_y_dict = parse_json_meta_data(
        reduction_input,
        "smearingPixelSizeY",
        1e-3,
        beam_center_run=True,
        background_run=True,
        empty_transmission_run=True,
        transmission_run=True,
        background_transmission=True,
        block_beam_run=True,
        dark_current_run=True,
    )

    # special loading case for sample to allow the slicing options
    logslice_data_dict = {}

    # Retrieve parameters for overwriting geometry related meta data
    swd_value_dict = parse_json_meta_data(
        reduction_input,
        "sampleToSi",
        1e-3,
        beam_center_run=True,
        background_run=True,
        empty_transmission_run=True,
        transmission_run=True,
        background_transmission=True,
        block_beam_run=True,
        dark_current_run=False,
    )
    # Sample to detector distance
    sdd_value_dict = parse_json_meta_data(
        reduction_input,
        "sampleDetectorDistance",
        1.0,
        beam_center_run=True,
        background_run=True,
        empty_transmission_run=True,
        transmission_run=True,
        background_transmission=True,
        block_beam_run=True,
        dark_current_run=False,
    )

    if timeslice or logslice:
        ws_name = f"{prefix}_{instrument_name}_{sample}_raw_histo_slice_group"
        if not registered_workspace(ws_name):
            filename = abspath(sample.strip(), instrument=instrument_name, ipts=ipts, directory=path)
            logger.notice(f"Loading filename {filename}")
            if timeslice:
                timesliceinterval = float(reduction_config["timeSliceInterval"])
                timesliceperiod = reduction_config["timeSlicePeriod"]
                logslicename = logsliceinterval = None
            elif logslice:
                timesliceinterval, timesliceperiod = None, None
                logslicename = reduction_config["logSliceName"]
                logsliceinterval = reduction_config["logSliceInterval"]
            biosans.load_and_split(
                filename,
                output_workspace=ws_name,
                time_interval=timesliceinterval,
                time_offset=reduction_config["timeSliceOffset"],
                time_period=timesliceperiod,
                log_name=logslicename,
                log_value_interval=logsliceinterval,
                sample_to_si_name=SAMPLE_SI_META_NAME,
                si_nominal_distance=SI_WINDOW_NOMINAL_DISTANCE_METER,
                sample_to_si_value=swd_value_dict[meta_data.SAMPLE],
                sample_detector_distance_value=sdd_value_dict[meta_data.SAMPLE],
                **load_params,
            )

            for _w in mtd[ws_name]:
                # Overwrite meta data
                set_meta_data(
                    str(_w),
                    wave_length=wave_length_dict[meta_data.SAMPLE],
                    wavelength_spread=wave_length_spread_dict[meta_data.SAMPLE],
                    sample_thickness=thickness,
                    sample_aperture_diameter=sample_aperture_diameter,
                    source_aperture_diameter=source_aperture_diameter,
                    smearing_pixel_size_x=smearing_pixel_size_x_dict[meta_data.SAMPLE],
                    smearing_pixel_size_y=smearing_pixel_size_y_dict[meta_data.SAMPLE],
                )
                # Transform X-axis to wave length with spread
                _w = transform_to_wavelength(_w)
                _w = set_init_uncertainties(_w)
                for btp_params in default_mask:
                    apply_mask(_w, **btp_params)

            if logslicename is not None:
                for n in range(mtd[ws_name].getNumberOfEntries()):
                    samplelogs = SampleLogs(mtd[ws_name].getItem(n))
                    logslice_data_dict[str(n)] = {
                        "data": list(samplelogs[logslicename].value),
                        "units": samplelogs[logslicename].units,
                        "name": logslicename,
                    }
    else:
        # Load raw without slicing
        ws_name = f"{prefix}_{instrument_name}_{sample}_raw_histo"
        if not registered_workspace(ws_name):
            # if sample is not an absolute path to nexus file, convert it to the absolute path
            filename = abspaths(sample, instrument=instrument_name, ipts=ipts, directory=path)
            # Pass load params to be used in LoadEventAsWorkspace2D
            load_params["XCenter"] = wave_length_dict[meta_data.SAMPLE]
            load_params["XWidth"] = wave_length_spread_dict[meta_data.SAMPLE]
            logger.notice(f"Loading filename {filename} from sample {sample}")
            biosans.load_events_and_histogram(
                filename,
                output_workspace=ws_name,
                sample_to_si_name=SAMPLE_SI_META_NAME,
                si_nominal_distance=SI_WINDOW_NOMINAL_DISTANCE_METER,
                sample_to_si_value=swd_value_dict[meta_data.SAMPLE],
                sample_detector_distance_value=sdd_value_dict[meta_data.SAMPLE],
                **load_params,
            )
            # Overwrite meta data
            set_meta_data(
                ws_name,
                wave_length=wave_length_dict[meta_data.SAMPLE],
                wavelength_spread=wave_length_spread_dict[meta_data.SAMPLE],
                sample_thickness=thickness,
                sample_aperture_diameter=sample_aperture_diameter,
                source_aperture_diameter=source_aperture_diameter,
                smearing_pixel_size_x=smearing_pixel_size_x_dict[meta_data.SAMPLE],
                smearing_pixel_size_y=smearing_pixel_size_y_dict[meta_data.SAMPLE],
            )

            # Apply mask
            for btp_params in default_mask:
                apply_mask(ws_name, **btp_params)

    reduction_input["logslice_data"] = logslice_data_dict

    # load all other files
    # for run_number in [center, bkgd, empty, sample_trans, bkgd_trans, blocked_beam]:
    for run_number, run_type in [
        (center, meta_data.BEAM_CENTER),
        (bkgd, meta_data.BACKGROUND),
        (empty, meta_data.EMPTY_TRANSMISSION),
        (sample_trans, meta_data.TRANSMISSION),
        (bkgd_trans, meta_data.TRANSMISSION_BACKGROUND),
        (blocked_beam, meta_data.BLOCK_BEAM),
    ]:
        if run_number:
            ws_name = f"{prefix}_{instrument_name}_{run_number}_raw_histo"
            if not registered_workspace(ws_name):
                filename = abspaths(run_number, instrument=instrument_name, ipts=ipts, directory=path)
                # Pass load params to be used in LoadEventAsWorkspace2D
                load_params["XCenter"] = wave_length_dict[run_type]
                load_params["XWidth"] = wave_length_spread_dict[run_type]
                logger.notice(f"Loading filename {filename}")
                biosans.load_events_and_histogram(
                    filename,
                    output_workspace=ws_name,
                    sample_to_si_name=SAMPLE_SI_META_NAME,
                    si_nominal_distance=SI_WINDOW_NOMINAL_DISTANCE_METER,
                    sample_to_si_value=swd_value_dict[run_type],
                    sample_detector_distance_value=sdd_value_dict[run_type],
                    **load_params,
                )
                # Set the wave length and wave length spread
                set_meta_data(
                    ws_name,
                    wave_length=wave_length_dict[run_type],
                    wavelength_spread=wave_length_spread_dict[run_type],
                    sample_thickness=None,
                    sample_aperture_diameter=None,
                    source_aperture_diameter=None,
                    smearing_pixel_size_x=smearing_pixel_size_x_dict[run_type],
                    smearing_pixel_size_y=smearing_pixel_size_y_dict[run_type],
                )
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)

    # load dark current
    # step 0: check if mid-range detector is used
    if reduction_input["has_midrange_detector"]:
        dark_current_file_midrange = reduction_config.get("darkMidrangeFileName", None)
    else:
        dark_current_file_midrange = "no_midrange_detector"
    # step 1: gather main and wing detector dark current
    dark_current_file_main = reduction_config.get("darkMainFileName", None)
    dark_current_file_wing = reduction_config.get("darkWingFileName", None)
    # step 2: logic
    # if and only if when none of the dark current files are None, we will load dark current
    # otherwise set all dark current to None
    if dark_current_file_main and dark_current_file_wing and dark_current_file_midrange:
        # main
        dark_current_main = dark_current_correction(
            dark_current_file_main,
            default_mask,
            instrument_name,
            ipts,
            load_params,
            path,
            prefix,
            wave_length_dict[meta_data.DARK_CURRENT],
            wave_length_spread_dict[meta_data.DARK_CURRENT],
            swd_value_dict[meta_data.DARK_CURRENT],
            sdd_value_dict[meta_data.DARK_CURRENT],
            smearing_pixel_size_x_dict[meta_data.DARK_CURRENT],
            smearing_pixel_size_y_dict[meta_data.DARK_CURRENT],
        )
        # wing
        dark_current_wing = dark_current_correction(
            dark_current_file_wing,
            default_mask,
            instrument_name,
            ipts,
            load_params,
            path,
            prefix,
            wave_length_dict[meta_data.DARK_CURRENT],
            wave_length_spread_dict[meta_data.DARK_CURRENT],
            swd_value_dict[meta_data.DARK_CURRENT],
            sdd_value_dict[meta_data.DARK_CURRENT],
            smearing_pixel_size_x_dict[meta_data.DARK_CURRENT],
            smearing_pixel_size_y_dict[meta_data.DARK_CURRENT],
        )
        # midrange
        if dark_current_file_midrange != "no_midrange_detector":
            dark_current_midrange = dark_current_correction(
                dark_current_file_midrange,
                default_mask,
                instrument_name,
                ipts,
                load_params,
                path,
                prefix,
                wave_length_dict[meta_data.DARK_CURRENT],
                wave_length_spread_dict[meta_data.DARK_CURRENT],
                swd_value_dict[meta_data.DARK_CURRENT],
                sdd_value_dict[meta_data.DARK_CURRENT],
                smearing_pixel_size_x_dict[meta_data.DARK_CURRENT],
                smearing_pixel_size_y_dict[meta_data.DARK_CURRENT],
            )
        else:
            dark_current_midrange = None
    else:
        dark_current_main = None
        dark_current_wing = None
        dark_current_midrange = None

    # load required processed_files
    # step 0: check if mid-range detector is used
    if reduction_input["has_midrange_detector"]:
        flood_file_midrange = reduction_config.get("sensitivityMidrangeFileName", None)
    else:
        flood_file_midrange = "no_midrange_detector"
    # step 1: gather main and wing detector processed_files
    flood_file_main = reduction_config.get("sensitivityMainFileName", None)
    flood_file_wing = reduction_config.get("sensitivityWingFileName", None)
    # step 2: logic
    # if and only if when none of the flood file are None, we will load teh sensitivity
    # workspace, otherwise set all sensitivity workspace to None
    if flood_file_main and flood_file_wing and flood_file_midrange:
        # main
        sensitivity_main_ws_name = f"{prefix}_main_sensitivity"
        if not registered_workspace(sensitivity_main_ws_name):
            logger.notice(f"Loading filename {flood_file_main}")
            load_sensitivity_workspace(flood_file_main, output_workspace=sensitivity_main_ws_name)
        # wing
        sensitivity_wing_ws_name = f"{prefix}_wing_sensitivity"
        if not registered_workspace(sensitivity_wing_ws_name):
            logger.notice(f"Loading filename {flood_file_wing}")
            load_sensitivity_workspace(flood_file_wing, output_workspace=sensitivity_wing_ws_name)
        # midrange
        if flood_file_midrange != "no_midrange_detector":
            sensitivity_midrange_ws_name = f"{prefix}_midrange_sensitivity"
            if not registered_workspace(sensitivity_midrange_ws_name):
                logger.notice(f"Loading filename {flood_file_midrange}")
                load_sensitivity_workspace(flood_file_midrange, output_workspace=sensitivity_midrange_ws_name)
        else:
            sensitivity_midrange_ws_name = None
    else:
        sensitivity_main_ws_name = None
        sensitivity_wing_ws_name = None
        sensitivity_midrange_ws_name = None

    mask_ws = None
    custom_mask_file = reduction_input["configuration"]["maskFileName"]
    if custom_mask_file is not None:
        mask_ws_name = f"{prefix}_mask"
        if not registered_workspace(mask_ws_name):
            logger.notice(f"Loading filename {custom_mask_file}")
            mask_ws = load_mask(custom_mask_file, output_workspace=mask_ws_name)
        else:
            mask_ws = mtd[mask_ws_name]

    if registered_workspace(f"{prefix}_{instrument_name}_{sample}_raw_histo_slice_group"):
        raw_sample_ws = mtd[f"{prefix}_{instrument_name}_{sample}_raw_histo_slice_group"]
        raw_sample_ws_list = [w for w in raw_sample_ws]
    else:
        raw_sample_ws_list = [mtd[f"{prefix}_{instrument_name}_{sample}_raw_histo"]]
    raw_bkgd_ws = mtd[f"{prefix}_{instrument_name}_{bkgd}_raw_histo"] if bkgd else None
    raw_blocked_ws = mtd[f"{prefix}_{instrument_name}_{blocked_beam}_raw_histo"] if blocked_beam else None
    raw_center_ws = mtd[f"{prefix}_{instrument_name}_{center}_raw_histo"]
    raw_empty_ws = mtd[f"{prefix}_{instrument_name}_{empty}_raw_histo"] if empty else None
    raw_sample_trans_ws = mtd[f"{prefix}_{instrument_name}_{sample_trans}_raw_histo"] if sample_trans else None
    raw_bkg_trans_ws = mtd[f"{prefix}_{instrument_name}_{bkgd_trans}_raw_histo"] if bkgd_trans else None
    sensitivity_main_ws = mtd[sensitivity_main_ws_name] if sensitivity_main_ws_name else None
    sensitivity_wing_ws = mtd[sensitivity_wing_ws_name] if sensitivity_wing_ws_name else None
    sensitivity_midrange_ws = mtd[sensitivity_midrange_ws_name] if sensitivity_midrange_ws_name else None

    if debug_output:
        for raw_sample in raw_sample_ws_list:
            plot_detector(
                input_workspace=str(raw_sample),
                filename=form_output_name(raw_sample),
                backend="mpl",
            )
        for ws in [
            raw_bkgd_ws,
            raw_center_ws,
            raw_empty_ws,
            raw_sample_trans_ws,
            raw_bkg_trans_ws,
            raw_blocked_ws,
            dark_current_main,
            dark_current_wing,
            dark_current_midrange,
        ]:
            if ws is not None:
                plot_detector(
                    input_workspace=str(ws),
                    filename=form_output_name(ws),
                    backend="mpl",
                    imshow_kwargs={"norm": LogNorm(vmin=1)},
                )

    return dict(
        sample=raw_sample_ws_list,
        background=raw_bkgd_ws,
        center=raw_center_ws,
        empty=raw_empty_ws,
        sample_transmission=raw_sample_trans_ws,
        background_transmission=raw_bkg_trans_ws,
        blocked_beam=raw_blocked_ws,
        dark_current_main=dark_current_main,
        dark_current_wing=dark_current_wing,
        dark_current_midrange=dark_current_midrange,
        sensitivity_main=sensitivity_main_ws,
        sensitivity_wing=sensitivity_wing_ws,
        sensitivity_midrange=sensitivity_midrange_ws,
        mask=mask_ws,
    )


def check_overlap_stitch_configuration(reduction_input: dict) -> None:
    """Validate the overlap stitch configuration depending on whether the midrange detector is present

    If the midrange detector is absent, there are two detectors and one overlap region, i.e. one Qmin/Qmax.
    If the midrange detector is present, there are three detectors and two overlap regions, i.e. two Qmin/Qmax.

    Parameters
    ----------
    reduction_input
        Reduction configuration parameters

    Raises
    -------
    ValueError
        If one of the overlap stitch boundary parameters has the wrong number of values.
    """
    reduction_config = reduction_input["configuration"]
    num_params_overlap_stitch = 1
    if reduction_input["has_midrange_detector"] and not reduction_config["overlapStitchIgnoreMidrange"]:
        num_params_overlap_stitch = 2
    params_overlap_stitch = ["overlapStitchQmin", "overlapStitchQmax"]
    if reduction_config["1DQbinType"] == "wedge":
        wedge_names = ["wedge1", "wedge2"]
        params_overlap_stitch = [wedge + boundary for wedge in wedge_names for boundary in params_overlap_stitch]
    for param in params_overlap_stitch:
        # these parameters are optional, only check values if not None
        if reduction_config[param] and len(reduction_config[param]) != num_params_overlap_stitch:
            has_midrange_str = "present" if reduction_input["has_midrange_detector"] else "not present"
            raise ValueError(
                f'Configuration "{param}" must have length {num_params_overlap_stitch} since midrange '
                f"detector {has_midrange_str}"
            )


def dark_current_correction(
    dark_current_file,
    default_mask,
    instrument_name,
    ipts,
    load_params,
    path,
    prefix,
    wavelength,
    wavelength_spread_user,
    user_sample_si_distance,
    user_sample_detector_distance,
    smearing_pixel_size_x,
    smearing_pixel_size_y,
):
    """Calculate the dark current correction

    Parameters
    ----------
    dark_current_file
    default_mask
    instrument_name
    ipts
    load_params
    path
    prefix
    wavelength
    wavelength_spread_user
    user_sample_si_distance: float, None
        user specified (overwriting) sample to Si-window distance in unit of meter
    user_sample_detector_distance: float, None
        user specified (overwriting) sample to detector distance in unit of meter
    smearing_pixel_size_x: float, None
        user specified smearing pixel size along X-direction
    smearing_pixel_size_y: float, None
        user specified smearing pixel size along Y-direction

    Returns
    -------
    ~mantid.api.MatrixWorkspace
       dark current correction workspace
    """
    run_number = extract_run_number(dark_current_file)
    ws_name = f"{prefix}_{instrument_name}_{run_number}_raw_histo"
    if not registered_workspace(ws_name):
        logger.notice(f"Loading filename {dark_current_file}")
        # identify to use exact given path to NeXus or use OnCat instead
        temp_name = abspath(dark_current_file, instrument=instrument_name, ipts=ipts, directory=path)
        if os.path.exists(temp_name):
            dark_current_file = temp_name
        # Pass load params to be used in LoadEventAsWorkspace2D
        load_params["XCenter"] = wavelength
        load_params["XWidth"] = wavelength_spread_user
        biosans.load_events_and_histogram(
            dark_current_file,
            output_workspace=ws_name,
            sample_to_si_name=SAMPLE_SI_META_NAME,
            si_nominal_distance=SI_WINDOW_NOMINAL_DISTANCE_METER,
            sample_to_si_value=user_sample_si_distance,
            sample_detector_distance_value=user_sample_detector_distance,
            **load_params,
        )
        # Set the wave length and wave length spread
        set_meta_data(
            ws_name,
            wave_length=wavelength,
            wavelength_spread=wavelength_spread_user,
            sample_thickness=None,
            sample_aperture_diameter=None,
            source_aperture_diameter=None,
            smearing_pixel_size_x=smearing_pixel_size_x,
            smearing_pixel_size_y=smearing_pixel_size_y,
        )
        for btp_params in default_mask:
            apply_mask(ws_name, **btp_params)
        dark_current = mtd[ws_name]
    else:
        dark_current = mtd[ws_name]
    return dark_current


def prepare_data_workspaces(
    data,
    center_x=None,
    center_y=None,
    center_y_wing=None,
    center_y_midrange=None,
    dark_current=None,
    flux_method=None,  # normalization (time/monitor)
    monitor_fail_switch=False,
    mask_ws=None,  # apply a custom mask from workspace
    mask_detector=None,  # main and/or wing and/or midrange
    mask_panel=None,  # mask back or front panel
    mask_btp=None,  # mask bank/tube/pixel
    solid_angle=True,
    sensitivity_workspace=None,
    output_workspace=None,
):
    r"""
    Given a " raw"data workspace, this function provides the following:

        - centers the detector
        - subtracts dark current
        - normalize by time or monitor
        - applies masks
        - corrects for solid angle
        - corrects for sensitivity

    All steps are optional. data, mask_ws, dark_current are either None
    or histogram workspaces. This function does not load any file.

    Parameters
    ----------
    data: ~mantid.dataobjects.Workspace2D
        raw workspace (histogram)
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    center_y_wing: float
        Move the center r of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    center_y_midrange: float
        Move the center of the detector vertically to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: ~mantid.dataobjects.Workspace2D
        histogram workspace containing the dark current measurement
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
    monitor_fail_switch: bool
        Resort to normalization by 'time' if 'monitor' was selected but no monitor counts are available
    mask_ws: ~mantid.dataobjects.Workspace2D
        Mask workspace
    mask_detector: str, list
        Name of one or more instrument components to mask: (e.g `detector1,wing_detector`)
    mask_panel: str
        Either 'front' or 'back' to mask whole front or back panel.
    mask_btp: dict
        Additional properties to Mantid's MaskBTP algorithm
    solid_angle: bool
        Apply the solid angle correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    output_workspace: str
        The output workspace name. If None will create data.name()+output_suffix

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = str(data)
        output_workspace = output_workspace.replace("_raw_histo", "") + "_processed_histo"

    mtd[str(data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    if center_x is not None and center_y is not None and center_y_wing is not None:
        # check whether:
        # (1) there is no midrange detector information provided or
        # (2) it is provided and the center_y_midrange should be there
        if (not has_midrange_detector(output_workspace)) or (
            has_midrange_detector(output_workspace) and center_y_midrange is not None
        ):
            biosans.center_detector(
                output_workspace,
                center_x=center_x,
                center_y=center_y,
                center_y_wing=center_y_wing,
                center_y_midrange=center_y_midrange,
            )

    # Dark current
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)

    # Normalization
    if str(flux_method).lower() == "monitor":
        try:
            normalize_by_monitor(output_workspace)
        except NoMonitorMetadataError as e:
            if monitor_fail_switch:
                logger.warning(f"{e}. Resorting to normalization by time")
                normalize_by_time(output_workspace)
            else:
                logger.warning(
                    '. Setting "normalizationResortToTime": True will cause the'
                    " reduction to normalize by time if monitor counts are not available"
                )
                raise
    elif str(flux_method).lower() == "time":
        normalize_by_time(output_workspace)
    else:
        logger.notice("No time or monitor normalization is carried out")

    # Mask either detector
    if mask_detector is not None:
        MaskDetectors(output_workspace, ComponentList=mask_detector)

    # Additional masks
    if mask_btp is None:
        mask_btp = dict()
    apply_mask(output_workspace, panel=mask_panel, mask=mask_ws, **mask_btp)

    # Solid angle
    if solid_angle:
        solid_angle_correction(output_workspace)

    # Sensitivity
    if sensitivity_workspace is not None:
        apply_sensitivity_correction(output_workspace, sensitivity_workspace=sensitivity_workspace)

    return mtd[output_workspace]


def process_single_configuration(
    sample_ws_raw: Union[str, Workspace2D],
    sample_trans_ws: Workspace2D = None,
    sample_trans_value: float = None,
    bkg_ws_raw: Workspace2D = None,
    bkg_trans_ws: Workspace2D = None,
    bkg_trans_value: float = None,
    blocked_ws_raw: Workspace2D = None,
    theta_dependent_transmission: bool = True,
    center_x: float = None,
    center_y: float = None,
    center_y_wing: float = None,
    center_y_midrange: float = None,
    dark_current: Workspace2D = None,
    flux_method: str = None,
    monitor_fail_switch: bool = False,
    mask_ws: Workspace2D = None,
    mask_detector: Union[str, List] = None,
    mask_panel: str = None,
    mask_btp: dict = None,
    solid_angle: bool = True,
    sensitivity_workspace: Workspace2D = None,
    thickness: float = 1.0,
    absolute_scale_method: str = "standard",
    absolute_scale: float = 1.0,
    keep_processed_workspaces: bool = True,
    output_workspace: str = None,
    output_prefix: str = "",
):
    r"""
    This function provides full data processing for a single experimental configuration,
    starting from workspaces (no data loading is happening inside this function)

    Steps applied to the sample and background workspaces:
    1. centering the detector
    2. dark current subtraction
    3. normalization by time or monitor
    4. masking
    5. solid angle correction
    6. sensitivity correction
    7. transmission correction

    Additional steps applied to the sample workspace:
    1. background subtraction
    2. thickness correction
    3. absolute scaling

    Parameters
    ----------
    sample_ws_raw
        raw data histogram workspace, or a string with the workspace name
    sample_trans_ws
        optional histogram workspace for sample transmission
    sample_trans_value
        optional value for sample transmission
    bkg_ws_raw
        optional raw histogram workspace for background
    bkg_trans_ws
        optional histogram workspace for background transmission
    bkg_trans_value
        optional value for background transmission
    blocked_ws_raw:
        optional histogram workspace for blocked beam
    theta_dependent_transmission: bool
        flag to apply angle dependent transmission
    center_x
        x center for the beam
    center_y
        y center for the beam
    center_y_wing
        y center for the wing detector
    center_y_midrange
        y center for the midrange detector
    dark_current
        dark current workspace
    flux_method
        normalization by 'time' or 'monitor'
    monitor_fail_switch
        resort to normalization by 'time' if 'monitor' was selected but no monitor counts are available
    mask_ws
        user-defined mask
    mask_detector
        Name of one or more instrument components to mask: (e.g `wing_detector,midrange_detector`)
    mask_panel
        mask front or back panel of tubes
    mask_btp
        optional bank, tube, pixel to mask. These are options passed as optional arguments to Mantid algorithm MaskBTP
    solid_angle
        flag to apply solid angle
    sensitivity_workspace
        workspace containing sensitivity
    thickness
        sample thickness (cm)
    absolute_scale_method
        method to do absolute scaling. One of ["standard", "direct_beam"]
    absolute_scale
        absolute scaling value for standard method
    keep_processed_workspaces
        flag to keep the processed blocked beam and background workspaces
    output_workspace
        output workspace name for the sample. If ``None``, the name will be ``output_suffix + "_sample"``
    output_prefix
        prefix for output workspaces (``output_suffix + "_sample"``, ``output_suffix + "_blocked"``, and
        ``output_suffix + "_background"``)

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = output_prefix + "_sample"

    # create a common configuration for prepare data
    prepare_data_conf = {
        "center_x": center_x,
        "center_y": center_y,
        "center_y_wing": center_y_wing,
        "center_y_midrange": center_y_midrange,
        "dark_current": dark_current,
        "flux_method": flux_method,
        "monitor_fail_switch": monitor_fail_switch,
        "mask_ws": mask_ws,
        "mask_detector": mask_detector,
        "mask_panel": mask_panel,
        "mask_btp": mask_btp,
        "solid_angle": solid_angle,
        "sensitivity_workspace": sensitivity_workspace,
    }

    # process the optional histogram workspace for blocked beam
    if blocked_ws_raw:
        blocked_ws_name = output_prefix + "_blocked"
        if not registered_workspace(blocked_ws_name):
            blocked_ws = prepare_data_workspaces(blocked_ws_raw, output_workspace=blocked_ws_name, **prepare_data_conf)
        else:
            blocked_ws = mtd[blocked_ws_name]

    # process sample
    sample_ws = prepare_data_workspaces(sample_ws_raw, output_workspace=output_workspace, **prepare_data_conf)
    if blocked_ws_raw:
        sample_ws = subtract_background(sample_ws, blocked_ws)
    # apply transmission to the sample
    transmission_dict = {}
    if sample_trans_ws or sample_trans_value:
        if sample_trans_value:
            transmission_dict = {"value": float(sample_trans_value), "error": ""}
        else:
            transmission_dict = {
                "value": sample_trans_ws.extractY(),
                "error": sample_trans_ws.extractE(),
            }
        sample_ws = apply_transmission_correction(
            sample_ws,
            trans_workspace=sample_trans_ws,
            trans_value=sample_trans_value,
            theta_dependent=theta_dependent_transmission,
            output_workspace=output_workspace,
        )

    # process background, if not already processed
    background_transmission_dict = {}
    if bkg_ws_raw:
        bkgd_ws_name = output_prefix + "_background"
        if not registered_workspace(bkgd_ws_name):
            bkgd_ws = prepare_data_workspaces(bkg_ws_raw, output_workspace=bkgd_ws_name, **prepare_data_conf)
            if blocked_ws_raw:
                bkgd_ws = subtract_background(bkgd_ws, blocked_ws)
            # apply transmission to bkgd
            if bkg_trans_ws or bkg_trans_value:
                if bkg_trans_value:
                    background_transmission_dict = {
                        "value": float(bkg_trans_value),
                        "error": "",
                    }
                else:
                    background_transmission_dict = {
                        "value": bkg_trans_ws.extractY(),
                        "error": bkg_trans_ws.extractE(),
                    }
                bkgd_ws = apply_transmission_correction(
                    bkgd_ws,
                    trans_workspace=bkg_trans_ws,
                    trans_value=bkg_trans_value,
                    theta_dependent=theta_dependent_transmission,
                    output_workspace=bkgd_ws_name,
                )
        else:
            bkgd_ws = mtd[bkgd_ws_name]
        # subtract background
        sample_ws = subtract_background(sample_ws, bkgd_ws, output_workspace=output_workspace)

        if not keep_processed_workspaces:
            bkgd_ws.delete()

    if blocked_ws_raw and not keep_processed_workspaces:
        blocked_ws.delete()

    # finalize with absolute scale and thickness
    sample_ws = normalize_by_thickness(sample_ws, thickness)

    # standard method assumes absolute scale from outside
    if absolute_scale_method == "direct_beam":
        raise NotImplementedError("This method is not yet implemented for BIOSANS")
    else:
        sample_ws *= absolute_scale

    return mtd[str(sample_ws)], {
        "sample": transmission_dict,
        "background": background_transmission_dict,
    }


def plot_reduction_output(
    reduction_output: list,
    reduction_input: dict,
    loglog: bool = True,
    imshow_kwargs: dict = {},
):
    """
    Generate reduction plot per slice per detector.

    Parameters
    ----------
    reduction_output:
        List of 1D and 2D I(Q) profile objects, for each detector panel
    reduction_input:.
        Dictionary of all possible options to properly configure the data reduction workflow.
    loglog:
        Whether to plot in loglog scale, by default True.
    imshow_kwargs:
        Keyword arguments to pass to imshow, by default {}.
    """
    reduction_config = reduction_input["configuration"]
    output_dir = reduction_config["outputDir"]
    outputFilename = reduction_input["outputFileName"]
    output_suffix = ""

    bin1d_type = reduction_config["1DQbinType"]

    imshow_kwargs = {} if imshow_kwargs is None else imshow_kwargs

    has_midrange_detector = reduction_input.get("has_midrange_detector", False)

    for i, out in enumerate(reduction_output):
        if len(reduction_output) > 1:
            output_suffix = f"_{i}"

        wedges = reduction_config["wedges"] if bin1d_type == "wedge" else None
        symmetric_wedges = reduction_config.get("symmetric_wedges", True)

        qmin_main = reduction_config["QminMain"]
        qmax_main = reduction_config["QmaxMain"]
        qmin_wing = reduction_config["QminWing"]
        qmax_wing = reduction_config["QmaxWing"]

        # special case for mid-range detector
        if has_midrange_detector:
            qmin_midrange = reduction_config["QminMidrange"]
            qmax_midrange = reduction_config["QmaxMidrange"]

        # NOTE: due to the api design in drtsans, pytest monkey patching does
        #       not work with the standard import, therefore we need to use
        #       full path import here to make unit test work.
        plot_IQazimuthal = drtsans.plots.api.plot_IQazimuthal
        plot_i1d = drtsans.plots.api.plot_i1d
        allow_overwrite = drtsans.path.allow_overwrite

        # main detector
        filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_main.png")
        plot_IQazimuthal(
            out.I2D_main,
            filename,
            backend="mpl",
            imshow_kwargs=imshow_kwargs,
            title="Main",
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            qmin=qmin_main,
            qmax=qmax_main,
        )
        plt.clf()

        # wing detector
        filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_wing.png")
        plot_IQazimuthal(
            out.I2D_wing,
            filename,
            backend="mpl",
            imshow_kwargs=imshow_kwargs,
            title="Wing",
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            qmin=qmin_wing,
            qmax=qmax_wing,
        )
        plt.clf()

        # special case for mid-range detector
        if has_midrange_detector:
            filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_midrange.png")
            plot_IQazimuthal(
                out.I2D_midrange,
                filename,
                backend="mpl",
                imshow_kwargs=imshow_kwargs,
                title="Midrange",
                wedges=wedges,
                symmetric_wedges=symmetric_wedges,
                qmin=qmin_midrange,
                qmax=qmax_midrange,
            )
            plt.clf()

        for j in range(len(out.I1D_main)):
            add_suffix = ""
            if len(out.I1D_main) > 1:
                add_suffix = f"_wedge_{j}"
            filename = os.path.join(output_dir, "1D", f"{outputFilename}{output_suffix}_1D{add_suffix}.png")
            if has_midrange_detector:
                plot_i1d(
                    [out.I1D_main[j], out.I1D_wing[j], out.I1D_midrange[j], out.I1D_combined[j]],
                    filename,
                    log_scale=loglog,
                    backend="mpl",
                    errorbar_kwargs={"label": "main,wing,midrange,both"},
                )
            else:
                plot_i1d(
                    [out.I1D_main[j], out.I1D_wing[j], out.I1D_combined[j]],
                    filename,
                    log_scale=loglog,
                    backend="mpl",
                    errorbar_kwargs={"label": "main,wing,both"},
                )
            plt.clf()

    plt.close()

    # allow overwrite
    allow_overwrite(os.path.join(output_dir, "1D"))
    allow_overwrite(os.path.join(output_dir, "2D"))


def reduce_single_configuration(
    loaded_ws: dict,
    reduction_input: dict,
    prefix: str = "",
    skip_nan: bool = True,
    debug_output: bool = False,
) -> list:
    """Do data reduction for the loaded workspaces according to the given reduction configuration

    Saves to file the resulting I(Q) profiles in ASCII format as well as a reduction log.

    Parameters
    ----------
    loaded_ws
        Dictionary with loaded workspaces (e.g. background, sensitivity, etc.) to use in the reduction
    reduction_input
        Reduction configuration parameters
    prefix
        String to prepend to output file names
    skip_nan
        If True, data points with value NaN will be ignored when writing I(Q) profiles to file.
    debug_output
        If True, saves 2D plots representative of the workspace at different stages of the reduction

    Returns
    -------
    list
        List of tuple of I(Q):s, one tuple per sample in the workspace. Each tuple holds I(Q):s with 1D and 2D binning
        for the individual detector panels, as well as a stitched 1D I(Q) profile.
    """
    reduction_config = reduction_input["configuration"]
    flux_method = reduction_config["normalization"]
    monitor_fail_switch = reduction_config["normalizationResortToTime"]
    transmission_radius = reduction_config["mmRadiusForTransmission"]
    solid_angle = reduction_config["useSolidAngleCorrection"]
    sample_trans_value = reduction_input["sample"]["transmission"]["value"]
    bkg_trans_value = reduction_input["background"]["transmission"]["value"]
    theta_dependent_transmission = reduction_config["useThetaDepTransCorrection"]
    mask_panel = None
    if reduction_config["useMaskBackTubes"] is True:
        mask_panel = "back"
    output_suffix = ""
    thickness = reduction_input["sample"]["thickness"]  # default thickness set in BIOSANS.json schema
    absolute_scale_method = reduction_config["absoluteScaleMethod"]
    absolute_scale = reduction_config["StandardAbsoluteScale"]
    time_slice = reduction_config["useTimeSlice"]
    time_slice_transmission = reduction_config["useTimeSlice"] and reduction_config["useTimeSliceTransmission"]
    # Transmission tolerance errors are only checked when slicing the transmission run,
    # because this is the only observed use-case of poor statistics for transmission calculations
    if time_slice_transmission:
        sample_trans_error_tol = reduction_input["sample"]["transmission"]["errorTolerance"]
        if sample_trans_error_tol is None:
            sample_trans_error_tol = DEFAULT_TRANSMISSION_ERROR_TOLERANCE
    else:
        sample_trans_error_tol = None
    output_dir = reduction_config["outputDir"]

    nxbins_main = reduction_config["numMainQxQyBins"]
    nybins_main = nxbins_main
    nxbins_wing = reduction_config["numWingQxQyBins"]
    nybins_wing = nxbins_wing
    nxbins_midrange = reduction_config["numMidrangeQxQyBins"]
    nybins_midrange = nxbins_midrange

    bin1d_type = reduction_config["1DQbinType"]
    log_binning = reduction_config["QbinType"] == "log"
    decade_on_center = reduction_config["useLogQBinsDecadeCenter"]  # default set in the schema

    nbins_main = reduction_config["numMainQBins"]
    nbins_main_per_decade = reduction_config["LogQBinsPerDecadeMain"]
    nbins_wing = reduction_config["numWingQBins"]
    nbins_wing_per_decade = reduction_config["LogQBinsPerDecadeWing"]
    nbins_midrange = reduction_config["numMidrangeQBins"]
    nbins_midrange_per_decade = reduction_config["LogQBinsPerDecadeMidrange"]

    outputFilename = reduction_input["outputFileName"]
    weighted_errors = reduction_config["useErrorWeighting"]
    qmin_main = reduction_config["QminMain"]
    qmax_main = reduction_config["QmaxMain"]
    qmin_wing = reduction_config["QminWing"]
    qmax_wing = reduction_config["QmaxWing"]
    qmin_midrange = reduction_config["QminMidrange"]
    qmax_midrange = reduction_config["QmaxMidrange"]
    wedge1_qmin_main = reduction_config["wedge1QminMain"]
    wedge1_qmax_main = reduction_config["wedge1QmaxMain"]
    wedge1_qmin_wing = reduction_config["wedge1QminWing"]
    wedge1_qmax_wing = reduction_config["wedge1QmaxWing"]
    wedge1_qmin_midrange = reduction_config["wedge1QminMidrange"]
    wedge1_qmax_midrange = reduction_config["wedge1QmaxMidrange"]
    wedge2_qmin_main = reduction_config["wedge2QminMain"]
    wedge2_qmax_main = reduction_config["wedge2QmaxMain"]
    wedge2_qmin_wing = reduction_config["wedge2QminWing"]
    wedge2_qmax_wing = reduction_config["wedge2QmaxWing"]
    wedge2_qmin_midrange = reduction_config["wedge2QminMidrange"]
    wedge2_qmax_midrange = reduction_config["wedge2QmaxMidrange"]
    annular_bin = reduction_config["AnnularAngleBin"]

    wedges_min = reduction_config["WedgeMinAngles"]
    wedges_max = reduction_config["WedgeMaxAngles"]
    wedges = None if wedges_min is None or wedges_max is None else list(zip(wedges_min, wedges_max))

    # automatically determine wedge binning if it wasn't explicitly set
    autoWedgeOpts = {}
    symmetric_wedges = True
    if bin1d_type == "wedge" and wedges_min is None:
        # the JSON validator "wedgesources" guarantees that the parameters to be collected are all non-empty
        autoWedgeOpts = {
            "q_min": reduction_config["autoWedgeQmin"],
            "q_delta": reduction_config["autoWedgeQdelta"],
            "q_max": reduction_config["autoWedgeQmax"],
            "azimuthal_delta": reduction_config["autoWedgeAzimuthalDelta"],
            "peak_width": reduction_config["autoWedgePeakWidth"],
            "background_width": reduction_config["autoWedgeBackgroundWidth"],
            "signal_to_noise_min": reduction_config["autoWedgeSignalToNoiseMin"],
            "auto_wedge_phi_min": reduction_config["autoWedgePhiMin"],
            "auto_wedge_phi_max": reduction_config["autoWedgePhiMax"],
            "auto_symmetric_wedges": reduction_config["autoSymmetricWedges"],
        }
        # auto-aniso returns all of the wedges
        symmetric_wedges = False

    fbc_options = biosans.fbc_options_json(reduction_input)
    xc, yc, yw, ym, fit_results = biosans.find_beam_center(loaded_ws.center, **fbc_options)
    logger.notice(f"Find beam center  = {xc}, {yc}, {yw}, {ym}")

    # does the run include the midrange detector? check the geometry of the first sample
    reduction_config["has_midrange_detector"] = has_midrange_detector(loaded_ws.sample[0])
    all_detectors_except_main = ["wing_detector"]
    all_detectors_except_wing = ["detector1"]
    if reduction_config["has_midrange_detector"]:
        all_detectors_except_main.append("midrange_detector")
        all_detectors_except_wing.append("midrange_detector")

    # empty beam transmission workspace
    if loaded_ws.empty is not None:
        empty_trans_ws_name = f"{prefix}_empty"
        empty_trans_ws = prepare_data_workspaces(
            loaded_ws.empty,
            flux_method=flux_method,
            monitor_fail_switch=monitor_fail_switch,
            mask_detector=all_detectors_except_main,
            center_x=xc,
            center_y=yc,
            center_y_wing=yw,
            center_y_midrange=ym,
            solid_angle=False,
            sensitivity_workspace=loaded_ws.sensitivity_main,
            output_workspace=empty_trans_ws_name,
        )
        if debug_output:
            plot_detector(
                empty_trans_ws,
                form_output_name(empty_trans_ws),
                backend="mpl",
                imshow_kwargs={"norm": LogNorm(vmin=1)},
            )
    else:
        empty_trans_ws = None

    # background transmission
    if loaded_ws.background_transmission is not None and empty_trans_ws is not None:
        bkgd_trans_ws_name = f"{prefix}_bkgd_trans"
        bkgd_trans_ws_processed = prepare_data_workspaces(
            loaded_ws.background_transmission,
            flux_method=flux_method,
            monitor_fail_switch=monitor_fail_switch,
            mask_detector=all_detectors_except_main,
            center_x=xc,
            center_y=yc,
            center_y_wing=yw,
            center_y_midrange=ym,
            solid_angle=False,
            sensitivity_workspace=loaded_ws.sensitivity_main,
            output_workspace=bkgd_trans_ws_name,
        )
        bkgd_trans_ws = calculate_transmission(
            bkgd_trans_ws_processed,
            empty_trans_ws,
            radius=transmission_radius,
            radius_unit="mm",
        )
        logger.notice(f"Background transmission = {bkgd_trans_ws.extractY()[0, 0]}")

        if debug_output:
            plot_detector(
                bkgd_trans_ws_processed,
                form_output_name(bkgd_trans_ws_processed),
                backend="mpl",
                imshow_kwargs={"norm": LogNorm(vmin=1)},
            )

    else:
        bkgd_trans_ws = None

    # sample transmission
    def _prepare_sample_transmission_ws(_sample_transmission):
        """
        inline function that prepare the sample transmission workspace for
        normalization usage
        """
        _ws_processed = prepare_data_workspaces(
            _sample_transmission,
            flux_method=flux_method,
            monitor_fail_switch=monitor_fail_switch,
            mask_detector=all_detectors_except_main,
            center_x=xc,
            center_y=yc,
            center_y_wing=yw,
            center_y_midrange=ym,
            solid_angle=False,
            sensitivity_workspace=loaded_ws.sensitivity_main,
            output_workspace=f"{prefix}_sample_trans",
        )

        return _ws_processed, calculate_transmission(
            _ws_processed,
            empty_trans_ws,
            radius=transmission_radius,
            radius_unit="mm",
            transmission_error_tolerance=sample_trans_error_tol,
        )

    sample_trans_ws = None
    if loaded_ws.sample_transmission is not None and empty_trans_ws is not None:
        sample_trans_ws_processed, sample_trans_ws = _prepare_sample_transmission_ws(loaded_ws.sample_transmission)
        logger.notice(f"Sample transmission = {sample_trans_ws.extractY()[0, 0]}")

        if debug_output:
            plot_detector(
                sample_trans_ws_processed,
                form_output_name(sample_trans_ws_processed),
                backend="mpl",
                imshow_kwargs={"norm": LogNorm(vmin=1)},
            )

    output = []
    detectordata = {}
    for i, raw_sample_ws in enumerate(loaded_ws.sample):
        sample_name = "_slice_{}".format(i + 1)
        if len(loaded_ws.sample) > 1:
            output_suffix = f"_{i}"

        if time_slice_transmission:
            try:
                _, sample_trans_ws = _prepare_sample_transmission_ws(raw_sample_ws)
            except (ZeroMonitorCountsError, TransmissionErrorToleranceError, TransmissionNanError) as e:
                logger.warning(f"Skipping slice {sample_name}: {e}.")
                continue

        try:
            processed_data_main, trans_main = process_single_configuration(
                raw_sample_ws,
                sample_trans_ws=sample_trans_ws,
                sample_trans_value=sample_trans_value,
                bkg_ws_raw=loaded_ws.background,
                bkg_trans_ws=bkgd_trans_ws,
                bkg_trans_value=bkg_trans_value,
                blocked_ws_raw=loaded_ws.blocked_beam,
                theta_dependent_transmission=theta_dependent_transmission,
                center_x=xc,
                center_y=yc,
                center_y_wing=yw,
                center_y_midrange=ym,
                dark_current=loaded_ws.dark_current_main,
                flux_method=flux_method,
                monitor_fail_switch=monitor_fail_switch,
                mask_detector=all_detectors_except_main,
                mask_ws=loaded_ws.mask,
                mask_panel=mask_panel,
                solid_angle=solid_angle,
                sensitivity_workspace=loaded_ws.sensitivity_main,
                output_workspace=f"processed_data_main_{i}",
                output_prefix=output_suffix,
                thickness=thickness,
                absolute_scale_method=absolute_scale_method,
                absolute_scale=absolute_scale,
                keep_processed_workspaces=False,
            )
            processed_data_wing, trans_wing = process_single_configuration(
                raw_sample_ws,
                sample_trans_ws=sample_trans_ws,
                sample_trans_value=sample_trans_value,
                bkg_ws_raw=loaded_ws.background,
                bkg_trans_ws=bkgd_trans_ws,
                bkg_trans_value=bkg_trans_value,
                blocked_ws_raw=loaded_ws.blocked_beam,
                theta_dependent_transmission=theta_dependent_transmission,
                center_x=xc,
                center_y=yc,
                center_y_wing=yw,
                center_y_midrange=ym,
                dark_current=loaded_ws.dark_current_wing,
                flux_method=flux_method,
                monitor_fail_switch=monitor_fail_switch,
                mask_detector=all_detectors_except_wing,
                mask_ws=loaded_ws.mask,
                mask_panel=mask_panel,
                solid_angle=solid_angle,
                sensitivity_workspace=loaded_ws.sensitivity_wing,
                output_workspace=f"processed_data_wing_{i}",
                output_prefix=output_suffix,
                thickness=thickness,
                absolute_scale_method=absolute_scale_method,
                absolute_scale=absolute_scale,
                keep_processed_workspaces=False,
            )
            if reduction_config["has_midrange_detector"]:
                processed_data_midrange, trans_midrange = process_single_configuration(
                    raw_sample_ws,
                    sample_trans_ws=sample_trans_ws,
                    sample_trans_value=sample_trans_value,
                    bkg_ws_raw=loaded_ws.background,
                    bkg_trans_ws=bkgd_trans_ws,
                    bkg_trans_value=bkg_trans_value,
                    blocked_ws_raw=loaded_ws.blocked_beam,
                    theta_dependent_transmission=theta_dependent_transmission,
                    center_x=xc,
                    center_y=yc,
                    center_y_wing=yw,
                    center_y_midrange=ym,
                    dark_current=loaded_ws.dark_current_midrange,
                    flux_method=flux_method,
                    monitor_fail_switch=monitor_fail_switch,
                    mask_detector=["detector1", "wing_detector"],
                    mask_ws=loaded_ws.mask,
                    mask_panel=mask_panel,
                    solid_angle=solid_angle,
                    sensitivity_workspace=loaded_ws.sensitivity_midrange,
                    output_workspace=f"processed_data_midrange_{i}",
                    output_prefix=output_suffix,
                    thickness=thickness,
                    absolute_scale_method=absolute_scale_method,
                    absolute_scale=absolute_scale,
                    keep_processed_workspaces=False,
                )
            else:
                processed_data_midrange, trans_midrange = (
                    None,
                    {
                        "sample": None,
                        "background": None,
                    },
                )
        except ZeroMonitorCountsError as e:
            if time_slice:
                logger.warning(f"Skipping slice {sample_name}: {e}.")
                continue
            else:
                raise

        if debug_output:
            from mantid.simpleapi import SaveNexusProcessed

            main_name = f"{form_output_name(processed_data_main).split('.')[0]}.nxs"
            wing_name = f"{form_output_name(processed_data_wing).split('.')[0]}.nxs"
            # remove history to write less data and speed up I/O
            if reduction_config["removeAlgorithmHistory"]:
                RemoveWorkspaceHistory(processed_data_main)
                RemoveWorkspaceHistory(processed_data_wing)
            SaveNexusProcessed(InputWorkspace=processed_data_main, Filename=main_name)
            SaveNexusProcessed(InputWorkspace=processed_data_wing, Filename=wing_name)
            plot_detector(
                processed_data_main,
                form_output_name(processed_data_main),
                backend="mpl",
            )  # imshow_kwargs={'norm': LogNorm(vmin=1)})
            # FIXME - using LogNorm option results in exception from matplotlib as some negative value detected
            plot_detector(
                processed_data_wing,
                form_output_name(processed_data_wing),
                backend="mpl",
            )  # , imshow_kwargs={'norm': LogNorm(vmin=1)})
            if reduction_config["has_midrange_detector"]:
                midrange_name = f"{form_output_name(processed_data_midrange).split('.')[0]}.nxs"
                # remove history to write less data and speed up I/O
                if reduction_config["removeAlgorithmHistory"]:
                    RemoveWorkspaceHistory(processed_data_midrange)
                SaveNexusProcessed(InputWorkspace=processed_data_midrange, Filename=midrange_name)
                plot_detector(
                    processed_data_midrange,
                    form_output_name(processed_data_midrange),
                    backend="mpl",
                )

        _loginfo = f"Transmission (main detector)@{str(raw_sample_ws)}:  {trans_main}\n"
        _loginfo += f"Transmission (midrange detector)@{str(raw_sample_ws)}: {trans_midrange}\n"
        _loginfo += f"Transmission (wing detector)@{str(raw_sample_ws)}: {trans_wing}"
        logger.notice(_loginfo)

        # binning
        subpixel_kwargs = dict()
        if reduction_config["useSubpixels"] is True:
            subpixel_kwargs = {
                "n_horizontal": reduction_config["subpixelsX"],
                "n_vertical": reduction_config["subpixelsY"],
            }
        iq1d_main_in = biosans.convert_to_q(processed_data_main, mode="scalar", **subpixel_kwargs)
        iq2d_main_in = biosans.convert_to_q(processed_data_main, mode="azimuthal", **subpixel_kwargs)
        if bool(autoWedgeOpts):  # determine wedges automatically from the main detector
            logger.notice(f"Auto wedge options: {autoWedgeOpts}")
            autoWedgeOpts["debug_dir"] = output_dir
            wedges = getWedgeSelection(iq2d_main_in, **autoWedgeOpts)
            logger.notice("found wedge angles:")
            peak_wedge, back_wedge = wedges
            logger.notice(f"    peak:      {peak_wedge}")
            logger.notice(f"    background: {back_wedge}")
            del peak_wedge, back_wedge
            logger.notice(f"wedges: {wedges}")

        # set the found wedge values to the reduction input, this will allow correct plotting
        reduction_config["wedges"] = wedges
        reduction_config["symmetric_wedges"] = symmetric_wedges

        iq2d_main_out, i1d_main_out = bin_all(
            iq2d_main_in,
            iq1d_main_in,
            nxbins_main,
            nybins_main,
            n1dbins=nbins_main,
            n1dbins_per_decade=nbins_main_per_decade,
            decade_on_center=decade_on_center,
            bin1d_type=bin1d_type,
            log_scale=log_binning,
            qmin=qmin_main,
            qmax=qmax_main,
            wedge1_qmin=wedge1_qmin_main,
            wedge1_qmax=wedge1_qmax_main,
            wedge2_qmin=wedge2_qmin_main,
            wedge2_qmax=wedge2_qmax_main,
            annular_angle_bin=annular_bin,
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            error_weighted=weighted_errors,
        )
        iq1d_wing_in = biosans.convert_to_q(processed_data_wing, mode="scalar")
        iq2d_wing_in = biosans.convert_to_q(processed_data_wing, mode="azimuthal")
        iq2d_wing_out, i1d_wing_out = bin_all(
            iq2d_wing_in,
            iq1d_wing_in,
            nxbins_wing,
            nybins_wing,
            n1dbins=nbins_wing,
            n1dbins_per_decade=nbins_wing_per_decade,
            decade_on_center=decade_on_center,
            bin1d_type=bin1d_type,
            log_scale=log_binning,
            qmin=qmin_wing,
            qmax=qmax_wing,
            wedge1_qmin=wedge1_qmin_wing,
            wedge1_qmax=wedge1_qmax_wing,
            wedge2_qmin=wedge2_qmin_wing,
            wedge2_qmax=wedge2_qmax_wing,
            annular_angle_bin=annular_bin,
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            error_weighted=weighted_errors,
        )
        if reduction_config["has_midrange_detector"]:
            iq1d_midrange_in = biosans.convert_to_q(processed_data_midrange, mode="scalar")
            iq2d_midrange_in = biosans.convert_to_q(processed_data_midrange, mode="azimuthal")
            iq2d_midrange_out, i1d_midrange_out = bin_all(
                iq2d_midrange_in,
                iq1d_midrange_in,
                nxbins_midrange,
                nybins_midrange,
                n1dbins=nbins_midrange,
                n1dbins_per_decade=nbins_midrange_per_decade,
                decade_on_center=decade_on_center,
                bin1d_type=bin1d_type,
                log_scale=log_binning,
                qmin=qmin_midrange,
                qmax=qmax_midrange,
                wedge1_qmin=wedge1_qmin_midrange,
                wedge1_qmax=wedge1_qmax_midrange,
                wedge2_qmin=wedge2_qmin_midrange,
                wedge2_qmax=wedge2_qmax_midrange,
                annular_angle_bin=annular_bin,
                wedges=wedges,
                symmetric_wedges=symmetric_wedges,
                error_weighted=weighted_errors,
            )
        else:
            iq2d_midrange_out, i1d_midrange_out = (IQazimuthal(intensity=[], error=[], qx=[], qy=[]), [])

        # create output directories
        create_output_dir(output_dir)

        # save ASCII and NXCANSAS files
        filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_main")
        save_ascii_binned_2D(f"{filename}.dat", "I(Qx,Qy)", iq2d_main_out)
        save_cansas_nx(iq2d_main_out.to_workspace(), f"{filename}.h5")
        filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_wing.dat")
        save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_wing_out)
        filename = os.path.join(output_dir, "2D", f"{outputFilename}{output_suffix}_2D_midrange.dat")
        save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_midrange_out)

        # intensity profiles in the order of stitching (lower to higher Q)
        if not reduction_config["has_midrange_detector"] or reduction_config["overlapStitchIgnoreMidrange"]:
            iq1d_profiles_in = [iq1d_main_in, iq1d_wing_in]
            i1d_profiles_out = [i1d_main_out, i1d_wing_out]
        else:
            iq1d_profiles_in = [iq1d_main_in, iq1d_midrange_in, iq1d_wing_in]
            i1d_profiles_out = [i1d_main_out, i1d_midrange_out, i1d_wing_out]
        iq1d_combined_out = stitch_binned_profiles(iq1d_profiles_in, i1d_profiles_out, reduction_config)

        save_i1d_all(
            i1d_main_out,
            i1d_midrange_out,
            i1d_wing_out,
            iq1d_combined_out,
            outputFilename,
            output_dir,
            output_suffix,
            skip_nan,
        )

        I_output = namedtuple(
            "I_output",
            ["I2D_main", "I2D_midrange", "I2D_wing", "I1D_main", "I1D_midrange", "I1D_wing", "I1D_combined"],
        )
        current_output = I_output(
            I2D_main=iq2d_main_out,
            I2D_midrange=iq2d_midrange_out,
            I2D_wing=iq2d_wing_out,
            I1D_main=i1d_main_out,
            I1D_midrange=i1d_midrange_out,
            I1D_wing=i1d_wing_out,
            I1D_combined=iq1d_combined_out,
        )
        output.append(current_output)

        detectordata[sample_name] = get_sample_detectordata(
            i1d_main_out,
            i1d_midrange_out,
            i1d_wing_out,
            iq1d_combined_out,
            iq2d_main_out,
            iq2d_midrange_out,
            iq2d_wing_out,
            reduction_config["has_midrange_detector"],
        )

    try:
        processed_data_main
    except NameError:
        raise NoDataProcessedError

    # save reduction log

    filename = os.path.join(
        reduction_config["outputDir"],
        outputFilename + f"_reduction_log{output_suffix}.hdf",
    )

    starttime = datetime.now().isoformat()
    reductionparams = {"data": copy.deepcopy(reduction_input)}
    specialparameters = {
        "beam_center": {"x": xc, "y": yc, "y_midrange": ym, "y_wing": yw},
        "fit results": fit_results,
        "sample_transmission": {
            "main": trans_main["sample"],
            "midrange": trans_midrange["sample"],
            "wing": trans_wing["sample"],
        },
        "background_transmission": {
            "main": trans_main["background"],
            "midrange": trans_midrange["background"],
            "wing": trans_wing["background"],
        },
    }

    samplelogs = {
        "main": SampleLogs(processed_data_main),
        "wing": SampleLogs(processed_data_wing),
    }
    if reduction_config["has_midrange_detector"]:
        samplelogs["midrange"] = SampleLogs(processed_data_midrange)
    logslice_data_dict = reduction_input["logslice_data"]

    savereductionlog(
        filename=filename,
        detectordata=detectordata,
        reductionparams=reductionparams,
        # pythonfile=pythonfile,
        starttime=starttime,
        specialparameters=specialparameters,
        logslicedata=logslice_data_dict,
        samplelogs=samplelogs,
    )

    # change permissions to all files to allow overwrite
    allow_overwrite(reduction_config["outputDir"])
    allow_overwrite(os.path.join(reduction_config["outputDir"], "1D"))
    allow_overwrite(os.path.join(reduction_config["outputDir"], "2D"))

    return output


def save_i1d_all(
    iq1d_main, iq1d_midrange, iq1d_wing, iq1d_combined, output_filename, output_dir, output_suffix, skip_nan
):
    """Save to file the intensity profiles for the individual detectors, as well as the combined (stitched) profile

    Parameters
    ----------
    iq1d_main: ~drtsans.dataobjects.IQmod
        Intensity profile for the main detector
    iq1d_midrange: ~drtsans.dataobjects.IQmod
        Intensity profile for the midrange detector
    iq1d_wing: ~drtsans.dataobjects.IQmod
        Intensity profile for the wing detector
    iq1d_combined: ~drtsans.dataobjects.IQmod
        Combined (stitched) intensity profile from different detectors
    output_filename: str
        The basename of the output files
    output_dir: str
        The directory to write files to
    output_suffix:str
        Suffix for output file names
    skip_nan: bool
        If true, any data point where intensity is NAN will not be written to file
    """
    for j in range(len(iq1d_main)):
        add_suffix = ""
        if len(iq1d_main) > 1:
            add_suffix = f"_wedge_{j}"
        main_1D_filename = os.path.join(
            output_dir,
            "1D",
            f"{output_filename}{output_suffix}_1D_main{add_suffix}.txt",
        )
        save_i1d(iq1d_main[j], main_1D_filename, skip_nan=skip_nan)

        wing_1D_filename = os.path.join(
            output_dir,
            "1D",
            f"{output_filename}{output_suffix}_1D_wing{add_suffix}.txt",
        )
        save_i1d(iq1d_wing[j], wing_1D_filename, skip_nan=skip_nan)

        if iq1d_midrange:
            midrange_1D_filename = os.path.join(
                output_dir,
                "1D",
                f"{output_filename}{output_suffix}_1D_midrange{add_suffix}.txt",
            )
            save_i1d(iq1d_midrange[j], midrange_1D_filename, skip_nan=skip_nan)

        combined_1D_filename = os.path.join(
            output_dir,
            "1D",
            f"{output_filename}{output_suffix}_1D_combined{add_suffix}.txt",
        )
        save_i1d(iq1d_combined[j], combined_1D_filename, skip_nan=skip_nan)


def get_sample_detectordata(
    iq1d_main, iq1d_midrange, iq1d_wing, iq1d_combined, iq2d_main, iq2d_midrange, iq2d_wing, has_midrange
):
    """Construct a dictionary of I(Q) output for the reduction logs

    Parameters
    ----------
    iq1d_main: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Intensity profile for the main detector
    iq1d_midrange: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Intensity profile for the midrange detector
    iq1d_wing: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Intensity profile for the wing detector
    iq1d_combined: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Combined (stitched) intensity profile from different detectors
    iq2d_main: ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy) for the main detector
    iq2d_midrange: ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy) for the midrange detector
    iq2d_wing: ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy) for the wing detector
    has_midrange: bool
        If True, the run includes the midrange detector
    """
    detector_data = {}
    # TODO: fix this bug - only one combined profile is saved in the reduction log, but for wedges there are two
    if iq1d_combined[-1].intensity.size > 0:
        detector_data["combined"] = {"i1d": [iq1d_combined[-1]]}

    # one 1D I(Q) per detector for scalar and annular binning and two 1D I(Q) per detector for wedge binning
    for index, (_iq1d_main, _iq1d_wing) in enumerate(zip(iq1d_main, iq1d_wing)):
        detector_data[f"main_{index}"] = {"i1d": [_iq1d_main]}
        detector_data[f"wing_{index}"] = {"i1d": [_iq1d_wing]}
    if has_midrange:
        for index, (_iq1d_midrange) in enumerate(iq1d_midrange):
            detector_data[f"midrange_{index}"] = {"i1d": [_iq1d_midrange]}

    # one 2D I(Q) per detector, assign to index 0
    index = 0
    detector_data[f"main_{index}"]["iqxqy"] = iq2d_main
    detector_data[f"wing_{index}"]["iqxqy"] = iq2d_wing
    if has_midrange:
        detector_data[f"midrange_{index}"]["iqxqy"] = iq2d_midrange

    return detector_data


def prepare_data(
    data,
    data_dir=None,
    scale_components=None,
    pixel_calibration=False,
    mask_detector=None,
    detector_offset=0,
    sample_offset=0,
    center_x=None,
    center_y=None,
    center_y_wing=None,
    center_y_midrange=None,
    dark_current=None,
    flux_method=None,
    monitor_fail_switch=False,
    mask=None,
    mask_panel=None,
    btp=dict(),
    solid_angle=True,
    sensitivity_file_path=None,
    sensitivity_workspace=None,
    wave_length=None,
    wavelength_spread=None,
    sample_aperture_diameter=None,
    sample_thickness=None,
    source_aperture_diameter=None,
    smearing_pixel_size_x=None,
    smearing_pixel_size_y=None,
    output_workspace=None,
    output_suffix="",
    **kwargs,
):
    r"""
    Load a BIOSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    data_dir: str, list
        Additional data search directories
    scale_components: Optional[dict]
        Dictionary of component names and scaling factors in the form of a three-element list,
        indicating rescaling of the pixels in the component along the X, Y, and Z axes.
        For instance, ``{"detector1": [1.0, 2.0, 1.0]}`` scales pixels along the Y-axis by a factor of 2,
        leaving the other pixel dimensions unchanged.
    pixel_calibration: bool
        Adjust pixel heights and widths according to barscan and tube-width calibrations.
    mask_detector: str, list
        Name of one or more instrument components to mask: (e.g `detector1,wing_detector`)
    detector_offset: float
        Additional translation of the detector along Z-axis, in mili-meters.
    sample_offset: float
        Additional translation of the sample along the Z-axis, in mili-meters.
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    center_y_wing: float
        Move the center of the wing detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    center_y_midrange: float
        Move the center of the detector vertically to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
    monitor_fail_switch: bool
        Resort to normalization by 'time' if 'monitor' was selected but no monitor counts are available
    panel: str
        Either 'front' or 'back' to mask a whole panel
    mask_panel: str
        Either 'front' or 'back' to mask whole front or back panel.
    mask: mask file path, MaskWorkspace, list
        Additional mask to be applied. If `list`, it is a list of
        detector ID's.
    btp: dict
        Additional properties to Mantid's MaskBTP algorithm
    solid_angle: bool
        Apply the solid angle correction
    sensitivity_file_path: str
        file containing previously calculated sensitivity correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    wave_length: float, None
        wave length in Angstrom
    wavelength_spread: float, None
        wave length spread in Angstrom
    sample_aperture_diameter: float, None
        sample aperture diameter in mm
    sample_thickness: None, float
        sample thickness in unit cm
    source_aperture_diameter: float, None
        source aperture size radius in unit mm
    smearing_pixel_size_x: float, None
        pixel size in x direction in unit as meter, only for Q-resolution calculation
    smearing_pixel_size_y: float, None
        pixel size in Y direction in unit as meter, only for Q-resolution calculation
    output_workspace: str
        Name of the output workspace. If not supplied, will be determined from the supplied value of ``data``.
    output_suffix: str
        If the ``output_workspace`` is not specified, this is appended to the automatically generated
        output workspace name.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    # Detector offset and sample offset are disabled
    if abs(detector_offset) > 1e-8 or abs(sample_offset) > 1e-8:
        raise RuntimeError("biosans.api.prepare_data does not work with detector_offset or sample_offset")

    # Load data and enforce to use nexus IDF
    if "enforce_use_nexus_idf" in kwargs:
        enforce_use_nexus_idf = kwargs["enforce_use_nexus_idf"]
    else:
        enforce_use_nexus_idf = False

    # Load event without moving detector and sample after loading NeXus and instrument
    ws = load_events(
        data,
        data_dir=data_dir,
        overwrite_instrument=True,
        output_workspace=output_workspace,
        output_suffix=output_suffix,
        scale_components=scale_components,
        pixel_calibration=pixel_calibration,
        detector_offset=0.0,
        sample_offset=0.0,
        LoadNexusInstrumentXML=enforce_use_nexus_idf,
    )

    # Reset the offset
    sample_offset, detector_offset = get_sample_detector_offset(
        ws, SAMPLE_SI_META_NAME, SI_WINDOW_NOMINAL_DISTANCE_METER
    )
    # Translate instrument with offsets
    move_instrument(ws, sample_offset, detector_offset)

    ws_name = str(ws)
    transform_to_wavelength(ws_name)
    set_init_uncertainties(ws_name)

    if center_x is not None and center_y is not None and center_y_wing is not None:
        # check whether:
        # (1) there is no midrange detector information provided or
        # (2) it is provided and the center_y_midrange should be there
        if (not has_midrange_detector(ws_name)) or (has_midrange_detector(ws_name) and center_y_midrange is not None):
            biosans.center_detector(
                ws_name,
                center_x=center_x,
                center_y=center_y,
                center_y_wing=center_y_wing,
                center_y_midrange=center_y_midrange,
            )

    # Mask either detector
    if mask_detector is not None:
        MaskDetectors(ws_name, ComponentList=mask_detector)

    # Dark current
    if dark_current is not None:
        if mtd.doesExist(str(dark_current)):
            dark_ws = mtd[str(dark_current)]
        else:
            # load dark current
            dark_ws = load_events(
                dark_current,
                overwrite_instrument=True,
                LoadNexusInstrumentXML=enforce_use_nexus_idf,
            )
            dark_ws = transform_to_wavelength(dark_ws)
            dark_ws = set_init_uncertainties(dark_ws)
        subtract_dark_current(ws_name, dark_ws)

    # Normalization
    if str(flux_method).lower() == "monitor":
        try:
            normalize_by_monitor(ws_name)
        except RuntimeError as e:
            if monitor_fail_switch:
                logger.warning(f"{e}. Resorting to normalization by time")
                normalize_by_time(ws_name)
            else:
                msg = (
                    '. Setting configuration "normalizationResortToTime": True will cause the'
                    " reduction to normalize by time if monitor counts are not available"
                )
                raise RuntimeError(str(e) + msg)
    elif str(flux_method).lower() == "time":
        normalize_by_time(ws_name)
    else:
        logger.notice("No time or monitor normalization is carried out")

    # Additional masks
    apply_mask(ws_name, panel=mask_panel, mask=mask, **btp)

    # Solid angle
    if solid_angle:
        if solid_angle is True:
            solid_angle_correction(ws_name)
        else:  # assume the solid_angle parameter is a workspace
            solid_angle_correction(ws_name, solid_angle_ws=solid_angle)
    # Sensitivity
    if sensitivity_file_path is not None or sensitivity_workspace is not None:
        drtsans.apply_sensitivity_correction(
            ws_name,
            sensitivity_filename=sensitivity_file_path,
            sensitivity_workspace=sensitivity_workspace,
        )

    # Overwrite meta data
    set_meta_data(
        ws_name,
        wave_length,
        wavelength_spread,
        sample_offset,
        sample_aperture_diameter,
        sample_thickness,
        source_aperture_diameter,
        smearing_pixel_size_x,
        smearing_pixel_size_y,
    )

    return mtd[ws_name]


def create_output_dir(output_dir):
    # output root
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # 1D and 2D
    for sub_dir in ["1D", "2D"]:
        n_d_dir = os.path.join(output_dir, sub_dir)
        if not os.path.exists(n_d_dir):
            os.mkdir(n_d_dir)


def form_output_name(workspace):
    workspace_name = str(workspace)
    file_name = workspace_name.split("/")[-1].split(".")[0]
    return f"{file_name}.png"


def file_has_midrange_detector(sample: str, instrument_name: str, ipts: str, directory: str) -> bool:
    """
    Check if the sample has midrange detector

    Parameters
    ----------
    sample
        Sample name
    instrument_name
        Instrument name
    ipts
        IPTS number
    directory
        Directory of the data

    Returns
    -------
    bool
        True if the sample has midrange detector
    """
    filename = abspaths(
        sample.strip(),
        instrument=instrument_name,
        ipts=ipts,
        directory=directory,
    ).split(",")[0]

    out_ws_name = mtd.unique_hidden_name()

    try:
        workspace = LoadEventNexus(Filename=filename, MetadataOnly=True, OutputWorkspace=out_ws_name)
    except RuntimeError:
        workspace = LoadNexusProcessed(Filename=filename, OutputWorkspace=out_ws_name)

    has_midrange = workspace.getInstrument().getComponentByName("midrange_detector") is not None

    # cleanup
    DeleteWorkspace(workspace)

    return has_midrange
