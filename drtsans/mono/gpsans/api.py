""" GPSANS API """
import os
import numpy as np
import ast
from collections import namedtuple

from mantid.simpleapi import mtd, MaskDetectors

from drtsans.path import registered_workspace
from drtsans.instruments import extract_run_number
from drtsans.settings import namedtuplefy
from drtsans.plots import plot_IQmod, plot_IQazimuthal
from drtsans.solid_angle import solid_angle_correction
from drtsans.beam_finder import center_detector, find_beam_center
from drtsans.mask_utils import apply_mask, load_mask
from drtsans.mono.load import load_events, transform_to_wavelength, load_events_and_histogram
from drtsans.mono.normalization import normalize_by_monitor, normalize_by_time
from drtsans.mono.dark_current import subtract_dark_current
from drtsans.sensitivity import apply_sensitivity_correction, load_sensitivity_workspace
from drtsans.transmission import apply_transmission_correction, calculate_transmission
from drtsans.thickness_normalization import normalize_by_thickness
from drtsans.mono.absolute_units import empty_beam_scaling
from drtsans.mono.gpsans.attenuation import attenuation_factor
from drtsans.mono.gpsans import convert_to_q
from drtsans import subtract_background
from drtsans.mono.meta_data import set_meta_data
from drtsans.iq import bin_all
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D


# Functions exposed to the general user (public) API
__all__ = ['prepare_data', 'prepare_data_workspaces', 'process_single_configuration',
           'load_all_files', 'plot_reduction_output', 'reduce_single_configuration']


@namedtuplefy
def load_all_files(reduction_input, prefix='', load_params=None):
    """
    load all required files at the beginning, and transform them to histograms
    """
    instrument_name = reduction_input["instrumentName"]
    ipts = reduction_input["iptsNumber"]
    sample = reduction_input["runNumber"]
    sample_trans = reduction_input["transmission"]["runNumber"]
    bkgd = reduction_input["background"]["runNumber"]
    bkgd_trans = reduction_input["background"]["transmission"]["runNumber"]
    empty = reduction_input["emptyTrans"]["runNumber"]
    center = reduction_input["beamCenter"]["runNumber"]
    if reduction_input["configuration"]["useBlockedBeam"]:
        blocked_beam = reduction_input["configuration"]["BlockBeamFileName"]
    else:
        blocked_beam = None

    # sample offsets, etc
    if load_params is None:
        load_params = {}
    wavelength = reduction_input["configuration"]["wavelength"]
    wavelengthSpread = reduction_input["configuration"]["wavelengthSpread"]
    if wavelength and wavelengthSpread:
        load_params["wavelengthinstrumentName"] = wavelength
        load_params["wavelengthSpread"] = wavelengthSpread

    if reduction_input["configuration"]["useDefaultMask"]:
        default_mask = [ast.literal_eval(mask_par) for mask_par in reduction_input["configuration"]["DefaultMask"]]
    else:
        default_mask = []
    for run_number in [center, sample, bkgd, empty, sample_trans, bkgd_trans, blocked_beam]:
        if run_number:
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                path = f"/HFIR/{instrument_name}/IPTS-{ipts}/nexus"
                filename = ','.join(f"{path}/{instrument_name}_{run.strip()}.nxs.h5" for run in run_number.split(','))
                print(f"Loading filename {filename}")
                load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)

    # do the same for dark current if exists
    dark_current = None
    if reduction_input["configuration"]["useDarkFileName"]:
        dark_current_file = reduction_input["configuration"]["darkFileName"]
        if dark_current_file:
            run_number = extract_run_number(dark_current_file)
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                print(f"Loading filename {dark_current_file}")
                dark_current = load_events_and_histogram(dark_current_file,
                                                         output_workspace=ws_name,
                                                         **load_params)
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)
            else:
                dark_current = mtd[ws_name]

    # load required processed_files
    sensitivity_ws_name = None
    if reduction_input["configuration"]["useSensitivityFileName"]:
        flood_file = reduction_input["configuration"]["sensitivityFileName"]
        if flood_file:
            sensitivity_ws_name = f'{prefix}_sensitivity'
            if not registered_workspace(sensitivity_ws_name):
                print(f"Loading filename {flood_file}")
                load_sensitivity_workspace(flood_file, output_workspace=sensitivity_ws_name)

    mask_ws = None
    if reduction_input["configuration"]["useMaskFileName"]:
        custom_mask_file = reduction_input["configuration"]["maskFileName"]
        if custom_mask_file:
            mask_ws_name = f'{prefix}_mask'
            if not registered_workspace(mask_ws_name):
                print(f"Loading filename {custom_mask_file}")
                mask_ws = load_mask(custom_mask_file, output_workspace=mask_ws_name)
            else:
                mask_ws = mtd[mask_ws_name]
    print('Done loading')

    raw_sample_ws = mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo']
    raw_bkgd_ws = mtd[f'{prefix}_{instrument_name}_{bkgd}_raw_histo'] if bkgd else None
    raw_blocked_ws = mtd[f'{prefix}_{instrument_name}_{blocked_beam}_raw_histo'] if blocked_beam else None
    raw_center_ws = mtd[f'{prefix}_{instrument_name}_{center}_raw_histo']
    raw_empty_ws = mtd[f'{prefix}_{instrument_name}_{empty}_raw_histo'] if empty else None
    raw_sample_trans_ws = mtd[f'{prefix}_{instrument_name}_{sample_trans}_raw_histo'] if sample_trans else None
    raw_bkg_trans_ws = mtd[f'{prefix}_{instrument_name}_{bkgd_trans}_raw_histo'] if bkgd_trans else None
    sensitivity_ws = mtd[sensitivity_ws_name] if sensitivity_ws_name else None

    return dict(sample=[raw_sample_ws],
                background=raw_bkgd_ws,
                center=raw_center_ws,
                empty=raw_empty_ws,
                sample_transmission=raw_sample_trans_ws,
                background_transmission=raw_bkg_trans_ws,
                blocked_beam=raw_blocked_ws,
                dark_current=dark_current,
                sensitivity=sensitivity_ws,
                mask=mask_ws)


def prepare_data(data,
                 mask_detector=None,
                 detector_offset=0, sample_offset=0,
                 center_x=None, center_y=None,
                 dark_current=None,
                 flux_method=None,
                 mask=None, mask_panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None, sensitivity_workspace=None,
                 wave_length=None, wavelength_spread=None,
                 sample_aperture_diameter=None, sample_thickness=None,
                 source_aperture_diameter=None,
                 pixel_size_x=None, pixel_size_y=None,
                 output_workspace=None, output_suffix='', **kwargs):
    r"""
    Load a GPSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    mask_detector: str
        Name of an instrument component to mask
    detector_offset: float
        Additional translation of the detector along Z-axis, in millimeters.
    sample_offset: float
        Additional translation of the sample along the Z-axis, in millimeters.
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
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
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter
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
    # GPSANS: detector offset is fixed to 0. Only detector sample distance is essential.
    #         So one offset is sufficient
    ws = load_events(data, overwrite_instrument=True, output_workspace=output_workspace, output_suffix=output_suffix,
                     detector_offset=0, sample_offset=sample_offset)

    ws_name = str(ws)
    transform_to_wavelength(ws_name)

    if center_x is not None and center_y is not None:
        center_detector(ws_name, center_x=center_x, center_y=center_y)

    # Dark current
    if dark_current is not None:
        if registered_workspace(str(dark_current)):
            dark_ws = mtd[str(dark_current)]
        else:
            dark_ws = load_events(dark_current, overwrite_instrument=True)
            dark_ws = transform_to_wavelength(dark_ws)
    else:
        dark_ws = None

    # load sensitivity
    if sensitivity_workspace is None and sensitivity_file_path:
        sensitivity_workspace = os.path.split(sensitivity_file_path)[-1]
        sensitivity_workspace = sensitivity_workspace.split('.')[0]
        sensitivity_workspace = load_sensitivity_workspace(sensitivity_file_path, sensitivity_workspace)

    # Mask either detector
    if mask_detector is not None:
        MaskDetectors(ws_name, ComponentList=mask_detector)

    # Overwrite meta data
    set_meta_data(ws_name, wave_length, wavelength_spread,
                  sample_offset,
                  sample_aperture_diameter, sample_thickness,
                  source_aperture_diameter,
                  pixel_size_x, pixel_size_y)

    return prepare_data_workspaces(ws_name,
                                   center_x=center_x,
                                   center_y=center_y,
                                   dark_current=dark_current,
                                   flux_method=flux_method,
                                   mask_ws=mask,
                                   mask_panel=mask_panel,
                                   mask_btp=btp,
                                   solid_angle=solid_angle,
                                   sensitivity_workspace=sensitivity_workspace,
                                   output_workspace=output_workspace,
                                   output_suffix=output_suffix)


def prepare_data_workspaces(data,
                            center_x=None, center_y=None,
                            dark_current=None,
                            flux_method=None,    # normalization (time/monitor)
                            mask_ws=None,        # apply a custom mask from workspace
                            mask_panel=None,     # mask back or front panel
                            mask_btp=None,       # mask bank/tube/pixel
                            solid_angle=True,
                            sensitivity_workspace=None,
                            output_workspace=None,
                            output_suffix='', **kwargs):
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
    dark_current: ~mantid.dataobjects.Workspace2D
        histogram workspace containing the dark current measurement
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
    mask_ws: ~mantid.dataobjects.Workspace2D
        Mask workspace
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
    output_suffix: str
        replace '_raw_histo' in the output workspace name.
        If empty, the default is '_processed_histo'

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = str(data)
        output_workspace.replace('_raw_histo', '') + '_processed_histo'

    mtd[str(data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    if center_x is not None and center_y is not None:
        center_detector(output_workspace, center_x=center_x, center_y=center_y)

    # Dark current
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)

    # Normalization
    if str(flux_method).lower() == 'monitor':
        normalize_by_monitor(output_workspace)
    elif str(flux_method).lower() == 'time':
        normalize_by_time(output_workspace)

    # Additional masks
    if mask_btp is None:
        mask_btp = dict()
    apply_mask(output_workspace, panel=mask_panel, mask=mask_ws, **mask_btp)

    # Solid angle
    if solid_angle:
        solid_angle_correction(output_workspace)

    # Sensitivity
    if sensitivity_workspace is not None:
        apply_sensitivity_correction(output_workspace,
                                     sensitivity_workspace=sensitivity_workspace)

    return mtd[output_workspace]


def process_single_configuration(sample_ws_raw,
                                 sample_trans_ws=None,
                                 sample_trans_value=None,
                                 bkg_ws_raw=None,
                                 bkg_trans_ws=None,
                                 bkg_trans_value=None,
                                 blocked_ws_raw=None,
                                 theta_deppendent_transmission=True,
                                 center_x=None,
                                 center_y=None,
                                 dark_current=None,
                                 flux_method=None,    # normalization (time/monitor)
                                 mask_ws=None,        # apply a custom mask from workspace
                                 mask_panel=None,     # mask back or front panel
                                 mask_btp=None,       # mask bank/tube/pixel
                                 solid_angle=True,
                                 sensitivity_workspace=None,
                                 output_workspace=None,
                                 output_suffix='',
                                 thickness=1.,
                                 absolute_scale_method='standard',
                                 empty_beam_ws=None,
                                 beam_radius=None,
                                 absolute_scale=1.,
                                 keep_processed_workspaces=True,
                                 **kwargs):
    r"""
    This function provides full data processing for a single experimental configuration,
    starting from workspaces (no data loading is happening inside this function)

    Parameters
    ----------
    sample_ws_raw: ~mantid.dataobjects.Workspace2D
        raw data histogram workspace
    sample_trans_ws: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for sample transmission
    sample_trans_value: float
        optional value for sample transmission
    bkg_ws_raw: ~mantid.dataobjects.Workspace2D
        optional raw histogram workspace for background
    bkg_trans_ws: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for background transmission
    bkg_trans_value: float
        optional value for background transmission
    blocked_ws_raw: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for blocked beam
    theta_deppendent_transmission: bool
        flag to apply angle dependent transmission
    center_x: float
        x center for the beam
    center_y: float
        y center for the beam
    dark_current: ~mantid.dataobjects.Workspace2D
        dark current workspace
    flux_method: str
        normalization by time or monitor
    mask_ws: ~mantid.dataobjects.Workspace2D
        user defined mask
    mask_panel: str
        mask fron or back panel
    mask_btp: dict
        optional bank, tube, pixel to mask
    solid_angle: bool
        flag to apply solid angle
    sensitivity_workspace: ~mantid.dataobjects.Workspace2D
        workspace containing sensitivity
    output_workspace: str
        output workspace name
    output_suffix:str
        suffix for output workspace
    thickness: float
        sample thickness (cm)
    absolute_scale_method: str
        method to do absolute scaling (standard or direct_beam)
    empty_beam_ws: ~mantid.dataobjects.Workspace2D
        empty beam workspace for absolute scaling
    beam_radius: float
        beam radius for absolute scaling
    absolute_scale: float
        absolute scaling value for standard method
    keep_processed_workspaces: bool
        flag to keep the processed blocked beam and background workspaces

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = output_suffix + '_sample'

    # create a common configuration for prepare data
    prepare_data_conf = {'center_x': center_x,
                         'center_y': center_y,
                         'dark_current': dark_current,
                         'flux_method': flux_method,
                         'mask_ws': mask_ws,
                         'mask_panel': mask_panel,
                         'mask_btp': mask_btp,
                         'solid_angle': solid_angle,
                         'sensitivity_workspace': sensitivity_workspace}

    # process blocked
    if blocked_ws_raw:
        blocked_ws_name = output_suffix + '_blocked'
        if not registered_workspace(blocked_ws_name):
            blocked_ws = prepare_data_workspaces(blocked_ws_raw,
                                                 output_workspace=blocked_ws_name,
                                                 **prepare_data_conf)
        else:
            blocked_ws = mtd[blocked_ws_name]

    # process sample
    sample_ws = prepare_data_workspaces(sample_ws_raw,
                                        output_workspace=output_workspace,
                                        **prepare_data_conf)
    if blocked_ws_raw:
        sample_ws = subtract_background(sample_ws, blocked_ws)
    # apply transmission to the sample
    if sample_trans_ws or sample_trans_value:
        sample_ws = apply_transmission_correction(sample_ws,
                                                  trans_workspace=sample_trans_ws,
                                                  trans_value=sample_trans_value,
                                                  theta_dependent=theta_deppendent_transmission,
                                                  output_workspace=output_workspace)

    # process background, if not already processed
    if bkg_ws_raw:
        bkgd_ws_name = output_suffix + '_background'
        if not registered_workspace(bkgd_ws_name):
            bkgd_ws = prepare_data_workspaces(bkg_ws_raw,
                                              output_workspace=bkgd_ws_name,
                                              **prepare_data_conf)
            if blocked_ws_raw:
                bkgd_ws = subtract_background(bkgd_ws, blocked_ws)
            # apply transmission to bkgd
            if bkg_trans_ws or bkg_trans_value:
                bkgd_ws = apply_transmission_correction(bkgd_ws,
                                                        trans_workspace=bkg_trans_ws,
                                                        trans_value=bkg_trans_value,
                                                        theta_dependent=theta_deppendent_transmission,
                                                        output_workspace=bkgd_ws_name)
        else:
            bkgd_ws = mtd[bkgd_ws_name]
        # subtract background
        sample_ws = subtract_background(sample_ws, bkgd_ws)

        if not keep_processed_workspaces:
            bkgd_ws.delete()

    if blocked_ws_raw and not keep_processed_workspaces:
        blocked_ws.delete()

    # finalize with absolute scale and thickness
    sample_ws = normalize_by_thickness(sample_ws, thickness)

    # standard method assumes absolute scale from outside
    if absolute_scale_method == 'direct_beam':
        try:
            empty = mtd[str(empty_beam_ws)]
        except KeyError:
            raise ValueError(f"Could not find empty beam {str(empty_beam_ws)}")

        ac, ace = attenuation_factor(empty)
        sample_ws = empty_beam_scaling(sample_ws,
                                       empty,
                                       beam_radius=beam_radius,
                                       unit='mm',
                                       attenuator_coefficient=ac,
                                       attenuator_error=ace,
                                       output_workspace=output_workspace)
    else:
        sample_ws *= absolute_scale

    return mtd[output_workspace]


def reduce_single_configuration(loaded_ws, reduction_input, prefix=''):
    flux_method = reduction_input["configuration"]["normalization"]
    try:
        transmission_radius = float(reduction_input["configuration"]["mmRadiusForTransmission"])
    except ValueError:
        transmission_radius = None
    solid_angle = reduction_input["configuration"]["useSolidAngleCorrection"]
    sample_trans_value = reduction_input["transmission"]["value"]
    bkg_trans_value = reduction_input["background"]["transmission"]["value"]
    theta_deppendent_transmission = reduction_input["configuration"]["useThetaDepTransCorrection"]
    mask_panel = None
    if reduction_input["configuration"]["useMaskBackTubes"]:
        mask_panel = 'back'
    output_suffix = ''
    try:
        thickness = float(reduction_input['thickness'])
    except ValueError:
        thickness = 1.
    absolute_scale_method = reduction_input["configuration"]["absoluteScaleMethod"]
    try:
        beam_radius = float(reduction_input["configuration"]["DBScalingBeamRadius"])
    except ValueError:
        beam_radius = None
    try:
        absolute_scale = float(reduction_input["configuration"]["StandardAbsoluteScale"])
    except ValueError:
        absolute_scale = 1.

    output_dir = reduction_input["configuration"]["outputDir"]

    nxbins_main = int(reduction_input["configuration"]["numQxQyBins"])
    nybins_main = int(nxbins_main)
    bin1d_type = reduction_input["configuration"]["1DQbinType"]
    log_binning = reduction_input["configuration"]["QbinType"] == 'log'
    even_decades = reduction_input["configuration"]["EvenDecades"]
    nbins_main = int(reduction_input["configuration"]["numQBins"])
    outputFilename = reduction_input["outputFilename"]
    weighted_errors = reduction_input["configuration"]["useErrorWeighting"]
    try:
        qmin = float(reduction_input["configuration"]["Qmin"])
    except ValueError:
        qmin = None
    qmax = None
    try:
        qmax = float(reduction_input["configuration"]["Qmax"])
    except ValueError:
        qmax = None
    try:
        annular_bin = float(reduction_input["configuration"]["AnnularAngleBin"])
    except ValueError:
        annular_bin = 1.
    wedges_min = np.fromstring(reduction_input["configuration"]["WedgeMinAngles"], sep=',')
    wedges_max = np.fromstring(reduction_input["configuration"]["WedgeMaxAngles"], sep=',')
    if len(wedges_min) != len(wedges_max):
        raise ValueError("The lengths of WedgeMinAngles and WedgeMaxAngles must be the same")
    wedges = list(zip(wedges_min, wedges_max))

    xc, yc = find_beam_center(loaded_ws.center)
    print("Center  =", xc, yc)

    # empty beam transmission workspace
    if loaded_ws.empty is not None:
        empty_trans_ws_name = f'{prefix}_empty'
        empty_trans_ws = prepare_data_workspaces(loaded_ws.empty,
                                                 flux_method=flux_method,
                                                 center_x=xc,
                                                 center_y=yc,
                                                 solid_angle=False,
                                                 sensitivity_ws=loaded_ws.sensitivity,
                                                 output_workspace=empty_trans_ws_name)
    else:
        empty_trans_ws = None

    # background transmission
    if loaded_ws.background_transmission:
        bkgd_trans_ws_name = f'{prefix}_bkgd_trans'
        bkgd_trans_ws_processed = prepare_data_workspaces(loaded_ws.background_transmission,
                                                          flux_method=flux_method,
                                                          center_x=xc,
                                                          center_y=yc,
                                                          solid_angle=False,
                                                          sensitivity_ws=loaded_ws.sensitivity,
                                                          output_workspace=bkgd_trans_ws_name)
        bkgd_trans_ws = calculate_transmission(bkgd_trans_ws_processed, empty_trans_ws,
                                               radius=transmission_radius, radius_unit="mm")
        print('Background transmission =', bkgd_trans_ws.extractY()[0, 0])
    else:
        bkgd_trans_ws = None

    # sample transmission
    if loaded_ws.sample_transmission:
        sample_trans_ws_name = f'{prefix}_sample_trans'
        sample_trans_ws_processed = prepare_data_workspaces(loaded_ws.sample_transmission,
                                                            flux_method=flux_method,
                                                            center_x=xc,
                                                            center_y=yc,
                                                            solid_angle=False,
                                                            sensitivity_ws=loaded_ws.sensitivity,
                                                            output_workspace=sample_trans_ws_name)
        sample_trans_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
                                                 radius=transmission_radius, radius_unit="mm")
        print('Sample transmission =', sample_trans_ws.extractY()[0, 0])
    else:
        sample_trans_ws = None

    output = []
    for i, raw_sample_ws in enumerate(loaded_ws.sample):
        if len(loaded_ws.sample) > 1:
            output_suffix = f'_{i}'
        processed_data_main = process_single_configuration(raw_sample_ws,
                                                           sample_trans_ws=sample_trans_ws,
                                                           sample_trans_value=sample_trans_value,
                                                           bkg_ws_raw=loaded_ws.background,
                                                           bkg_trans_ws=bkgd_trans_ws,
                                                           bkg_trans_value=bkg_trans_value,
                                                           blocked_ws_raw=loaded_ws.blocked_beam,
                                                           theta_deppendent_transmission=theta_deppendent_transmission,
                                                           center_x=xc, center_y=yc,
                                                           dark_current=loaded_ws.dark_current,
                                                           flux_method=flux_method,
                                                           mask_ws=loaded_ws.mask,
                                                           mask_panel=mask_panel,
                                                           solid_angle=solid_angle,
                                                           sensitivity_workspace=loaded_ws.sensitivity,
                                                           output_workspace=f'processed_data_main',
                                                           output_suffix=output_suffix,
                                                           thickness=thickness,
                                                           absolute_scale_method=absolute_scale_method,
                                                           empty_beam_ws=empty_trans_ws,
                                                           beam_radius=beam_radius,
                                                           absolute_scale=absolute_scale,
                                                           keep_processed_workspaces=False)
        # binning
        iq2d_main_in = convert_to_q(processed_data_main, mode='azimuthal')
        iq1d_main_in = convert_to_q(processed_data_main, mode='scalar')
        iq2d_main_out, iq1d_main_out = bin_all(iq2d_main_in, iq1d_main_in, nxbins_main, nybins_main, nbins_main,
                                               bin1d_type=bin1d_type, log_scale=log_binning,
                                               even_decade=even_decades, qmin=qmin, qmax=qmax,
                                               annular_angle_bin=annular_bin, wedges=wedges,
                                               error_weighted=weighted_errors)

        # save ASCII files
        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D.txt')
        save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_main_out)

        for j in range(len(iq1d_main_out)):
            add_suffix = ""
            if len(iq1d_main_out) > 1:
                add_suffix = f'_wedge_{j}'
            ascii_1D_filename = os.path.join(output_dir, '1D',
                                             f'{outputFilename}{output_suffix}_1D{add_suffix}.txt')
            save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq1d_main_out[j])

        IofQ_output = namedtuple('IofQ_output', ['I2D_main', 'I1D_main'])
        current_output = IofQ_output(I2D_main=iq2d_main_out,
                                     I1D_main=iq1d_main_out)
        output.append(current_output)
    return output


def plot_reduction_output(reduction_output, reduction_input, imshow_kwargs=None):
    output_dir = reduction_input["configuration"]["outputDir"]
    outputFilename = reduction_input["outputFilename"]
    output_suffix = ''
    if imshow_kwargs is None:
        imshow_kwargs = {}
    for i, out in enumerate(reduction_output):
        if len(reduction_output) > 1:
            output_suffix = f'_{i}'
        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D.png')
        plot_IQazimuthal(out.I2D_main, filename, backend='mpl', imshow_kwargs=imshow_kwargs, title='Main')
        for j in range(len(out.I1D_main)):
            add_suffix = ""
            if len(out.I1D_main) > 1:
                add_suffix = f'_wedge_{j}'
            filename = os.path.join(output_dir, '1D', f'{outputFilename}{output_suffix}_1D{add_suffix}.png')
            plot_IQmod([out.I1D_main[j]], filename, loglog=True,
                       backend='mpl', errorbar_kwargs={'label': 'main'})
