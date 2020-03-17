""" Top-level API for EQSANS """
from datetime import datetime
import copy
import numpy as np
import os

from collections import namedtuple
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

from mantid.simpleapi import mtd, logger, SaveAscii, SaveNexus  # noqa E402
# Import rolled up to complete a single top-level API
import drtsans  # noqa E402
from drtsans import (apply_sensitivity_correction, getWedgeSelection, load_sensitivity_workspace, solid_angle_correction)  # noqa E402
from drtsans import subtract_background  # noqa E402
from drtsans.settings import namedtuplefy  # noqa E402
from drtsans.process_uncertainties import set_init_uncertainties  # noqa E402
from drtsans.save_ascii import save_ascii_1D, save_xml_1D, save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.save_2d import save_nist_dat, save_nexus  # noqa E402
from drtsans.transmission import apply_transmission_correction  # noqa E402
from drtsans.tof.eqsans.transmission import calculate_transmission  # noqa E402
from drtsans.thickness_normalization import normalize_by_thickness  # noqa E402
from drtsans.beam_finder import find_beam_center  # noqa E402
from drtsans.instruments import extract_run_number  # noqa E402
from drtsans.path import registered_workspace  # noqa E402
from drtsans.tof.eqsans.load import load_events, load_events_and_histogram, load_and_split  # noqa E402
from drtsans.tof.eqsans.dark_current import subtract_dark_current  # noqa E402
from drtsans.tof.eqsans.cfg import load_config  # noqa E402
from drtsans.samplelogs import SampleLogs  # noqa E402
from drtsans.mask_utils import apply_mask, load_mask  # noqa E402
from drtsans.tof.eqsans.normalization import normalize_by_flux  # noqa E402
from drtsans.tof.eqsans.meta_data import set_meta_data  # noqa E402
from drtsans.tof.eqsans.momentum_transfer import convert_to_q, split_by_frame  # noqa E402
from drtsans.plots import plot_IQmod, plot_IQazimuthal  # noqa E402
from drtsans.iq import bin_all  # noqa E402
from drtsans.dataobjects import save_iqmod  # noqa E402

__all__ = ['apply_solid_angle_correction', 'subtract_background',
           'prepare_data', 'save_ascii_1D', 'save_xml_1D',
           'save_nist_dat', 'save_nexus', 'set_init_uncertainties',
           'load_all_files', 'prepare_data_workspaces',
           'process_single_configuration', 'reduce_single_configuration',
           'plot_reduction_output']


def _get_configuration_file_parameters(sample_run):
    try:
        configuration_file_parameters = load_config(source=sample_run)
    except RuntimeError as e:
        logger.error(e)
        logger.warning('Not using previous configuration')
        configuration_file_parameters = {}
    return configuration_file_parameters


@namedtuplefy
def load_all_files(reduction_input, prefix='', load_params=None):
    r"""
    overwrites metadata for sample workspace
    """
    instrument_name = reduction_input["instrumentName"]
    ipts = reduction_input["iptsNumber"]
    path = f'/SNS/{instrument_name}/IPTS-{ipts}/nexus'
    sample = reduction_input["runNumber"]
    sample_trans = reduction_input["transmission"]["runNumber"]
    bkgd = reduction_input["background"]["runNumber"]
    bkgd_trans = reduction_input["background"]["transmission"]["runNumber"]
    empty = reduction_input["empty"]["runNumber"]
    center = reduction_input["beamCenter"]["runNumber"]

    default_mask = None
    if reduction_input["configuration"]["useDefaultMask"]:
        configuration_file_parameters = _get_configuration_file_parameters(sample.split(',')[0].strip())
        default_mask = configuration_file_parameters['combined mask']

    # find the center first
    center = reduction_input["beamCenter"]["runNumber"]
    if center != "":
        center_ws_name = f'{prefix}_{instrument_name}_{center}_raw_events'
        if not registered_workspace(center_ws_name):
            center_filename = f'{path}/{instrument_name}_{center}.nxs.h5'
            load_events(center_filename, output_workspace=center_ws_name)
            if reduction_input["configuration"]["useDefaultMask"]:
                apply_mask(center_ws_name, mask=default_mask)
        center_x, center_y = find_beam_center(center_ws_name)
        logger.notice(f"calculated center ({center_x}, {center_y})")
        print(f"calculated center ({center_x}, {center_y})")
    else:
        center_x = 0.025239
        center_y = 0.0170801
        logger.notice(f"use default center ({center_x}, {center_y})")

    if load_params is None:
        load_params = dict(center_x=center_x, center_y=center_y, keep_events=False)
    if reduction_input["configuration"]["useDetectorOffset"] and reduction_input["configuration"]['detectorOffset']:
        load_params['detector_offset'] = float(reduction_input["configuration"]['detectorOffset'])
    if reduction_input["configuration"]["useSampleOffset"] and reduction_input["configuration"]['sampleOffset']:
        load_params['sample_offset'] = float(reduction_input["configuration"]['sampleOffset'])
    if reduction_input["configuration"]["useTOFcuts"] and reduction_input["configuration"]["TOFmin"].strip() != "":
        load_params['low_tof_clip'] = float(reduction_input["configuration"]['TOFmin'])
    else:
        load_params['low_tof_clip'] = 500.
    if reduction_input["configuration"]["useTOFcuts"] and reduction_input["configuration"]["TOFmax"].strip() != "":
        load_params['high_tof_clip'] = float(reduction_input["configuration"]['TOFmax'])
    else:
        load_params['high_tof_clip'] = 2000.
    if reduction_input["configuration"]["wavelengthStep"].strip() != "":
        load_params['bin_width'] = float(reduction_input["configuration"]['wavelengthStep'])
    load_params['monitors'] = reduction_input["configuration"]["normalization"] == "Monitor"

    # for now, there are some big issues with the monitor on EQSANS
    # since we cannot guarantee that it will work, this option is disabled for now
    if load_params['monitors']:
        raise RuntimeError('Normalization by monitor option will be enabled in a later drt-sans release')

    path = f"/SNS/{instrument_name}/IPTS-{ipts}/nexus"

    # check for time/log slicing
    timeslice = reduction_input["configuration"].get("timeslice")
    logslice = reduction_input["configuration"].get("logslice")

    if timeslice and logslice:
        raise ValueError("Can't do both time slicing and log slicing")

    if timeslice or logslice:
        if len(sample.split(',')) > 1:
            raise ValueError("Can't do slicing on summed data sets")

    # special loading case for sample to allow the slicing options
    if timeslice or logslice:
        ws_name = f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group'
        if not registered_workspace(ws_name):
            filename = f"{path}/{instrument_name}_{sample.strip()}.nxs.h5"
            print(f"Loading filename {filename}")
            if timeslice:
                timesliceinterval = float(reduction_input["configuration"]["timesliceinterval"])
                logslicename = None
                logsliceinterval = None
            elif logslice:
                timesliceinterval = None
                logslicename = reduction_input["configuration"]["logslicename"]
                logsliceinterval = float(reduction_input["configuration"]["logsliceinterval"])
            load_and_split(filename, output_workspace=ws_name,
                           time_interval=timesliceinterval,
                           log_name=logslicename, log_value_interval=logsliceinterval,
                           **load_params)
            for _w in mtd[ws_name]:
                if default_mask:
                    apply_mask(_w, mask=default_mask)
    else:
        ws_name = f'{prefix}_{instrument_name}_{sample}_raw_histo'
        if not registered_workspace(ws_name):
            filename = ','.join(f"{path}/{instrument_name}_{run.strip()}.nxs.h5" for run in sample.split(','))
            print(f"Loading filename {filename}")
            load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
            if default_mask:
                apply_mask(ws_name, mask=default_mask)

    # load all other files
    for run_number in [bkgd, empty, sample_trans, bkgd_trans]:
        if run_number:
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                filename = ','.join(f"{path}/{instrument_name}_{run.strip()}.nxs.h5" for run in run_number.split(','))
                print(f"Loading filename {filename}")
                load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
                if default_mask:
                    apply_mask(ws_name, mask=default_mask)

    dark_current_ws = None
    dark_current_mon_ws = None
    if reduction_input["configuration"]["useDarkFileName"]:
        dark_current_file = reduction_input["configuration"]["darkFileName"]
        if dark_current_file:
            run_number = extract_run_number(dark_current_file)
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                print(f"Loading filename {dark_current_file}")
                dark_current_ws, _ = load_events_and_histogram(dark_current_file,
                                                               output_workspace=ws_name,
                                                               **load_params)
                if default_mask:
                    apply_mask(ws_name, mask=default_mask)
            else:
                dark_current_ws = mtd[ws_name]

    if registered_workspace(f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group'):
        sample_ws = mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group']
        sample_ws_list = [w for w in sample_ws]
    else:
        sample_ws_list = [mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo']]
    background_ws = mtd[f'{prefix}_{instrument_name}_{bkgd}_raw_histo'] if bkgd else None
    empty_ws = mtd[f'{prefix}_{instrument_name}_{empty}_raw_histo'] if empty else None
    sample_transmission_ws = mtd[f'{prefix}_{instrument_name}_{sample_trans}_raw_histo'] if sample_trans else None
    background_transmission_ws = mtd[f'{prefix}_{instrument_name}_{bkgd_trans}_raw_histo'] if bkgd_trans else None
    if load_params['monitors']:
        sample_mon_ws = mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo_monitors']
        background_mon_ws = mtd[f'{prefix}_{instrument_name}_{bkgd}_raw_histo_monitors'] if bkgd else None
        empty_mon_ws = mtd[f'{prefix}_{instrument_name}_{empty}_raw_histo_monitors'] if empty else None
        sample_transmission_mon_ws = mtd[f'{prefix}_{instrument_name}_{sample_trans}' +
                                         f'_raw_histo_monitors'] if sample_trans else None
        background_transmission_mon_ws = mtd[f'{prefix}_{instrument_name}_{bkgd_trans}' +
                                             f'_raw_histo_monitors'] if bkgd_trans else None
    else:
        sample_mon_ws = None
        background_mon_ws = None
        empty_mon_ws = None
        sample_transmission_mon_ws = None
        background_transmission_mon_ws = None

    # load required processed_files
    sensitivity_ws = None
    if reduction_input["configuration"]["useSensitivityFileName"]:
        flood_file = reduction_input["configuration"]["sensitivityFileName"]
        if flood_file:
            sensitivity_ws_name = f'{prefix}_sensitivity'
            if not registered_workspace(sensitivity_ws_name):
                print(f"Loading filename {flood_file}")
                load_sensitivity_workspace(flood_file, output_workspace=sensitivity_ws_name)
            sensitivity_ws = mtd[sensitivity_ws_name]

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

    # TODO load these files only once
    # beam_flux_ws = None
    # monitor_flux_ratio_ws = None

    print('Done loading')
    try:
        sample_aperture_diameter = float(reduction_input["configuration"]["sampleApertureSize"])
    except (KeyError, ValueError):
        sample_aperture_diameter = None
    sample_thickness = float(reduction_input["thickness"])
    try:
        pixel_size_x = float(reduction_input["configuration"]["pixel_size_x"])
    except (KeyError, ValueError):
        pixel_size_x = None
    try:
        pixel_size_y = float(reduction_input["configuration"]["pixel_size_y"])
    except (KeyError, ValueError):
        pixel_size_y = None

    for ws in sample_ws_list:
        set_meta_data(ws, wave_length=None, wavelength_spread=None,
                      sample_offset=load_params['sample_offset'],
                      sample_aperture_diameter=sample_aperture_diameter,
                      sample_thickness=sample_thickness,
                      source_aperture_diameter=None,
                      pixel_size_x=pixel_size_x,
                      pixel_size_y=pixel_size_y)
    ws_mon_pair = namedtuple('ws_mon_pair', ['data', 'monitor'])
    return dict(sample=[ws_mon_pair(data=ws, monitor=sample_mon_ws) for ws in sample_ws_list],
                background=ws_mon_pair(data=background_ws, monitor=background_mon_ws),
                empty=ws_mon_pair(data=empty_ws, monitor=empty_mon_ws),
                sample_transmission=ws_mon_pair(data=sample_transmission_ws, monitor=sample_transmission_mon_ws),
                background_transmission=ws_mon_pair(data=background_transmission_ws,
                                                    monitor=background_transmission_mon_ws),
                dark_current=ws_mon_pair(data=dark_current_ws, monitor=dark_current_mon_ws),
                sensitivity=sensitivity_ws, mask=mask_ws)


def prepare_data_workspaces(data,
                            dark_current=None,
                            flux_method=None,    # normalization (proton charge/time/monitor)
                            flux=None,           # additional file for normalization
                            mask_ws=None,        # apply a custom mask from workspace
                            mask_panel=None,     # mask back or front panel
                            mask_btp=None,       # mask bank/tube/pixel
                            solid_angle=True,
                            sensitivity_workspace=None,
                            output_workspace=None,
                            output_suffix='', **kwargs):

    r"""
    Given a " raw"data workspace, this function provides the following:

        - subtracts dark current
        - normalize by time or monitor
        - applies masks
        - corrects for solid angle
        - corrects for sensitivity

    All steps are optional. data, mask_ws, dark_current are either None
    or histogram workspaces. This function does not load any file.

    Parameters
    ----------
    data: namedtuple
        (~mantid.dataobjects.Workspace2D, ~mantid.dataobjects.Workspace2D)
        raw workspace (histogram) for data and monitor
    dark_current: ~mantid.dataobjects.Workspace2D
        histogram workspace containing the dark current measurement
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
    flux: str
        if ``flux_method`` is proton charge, then path to file containing the
        wavelength distribution of the neutron flux. If ``flux method`` is
        monitor, then path to file containing the flux-to-monitor ratios.
        if ``flux_method`` is time, then pass one log entry name such
        as ``duration`` or leave it as :py:obj:`None` for automatic log search.
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
        output_workspace = str(data.data)
        output_workspace.replace('_raw_histo', '') + '_processed_histo'

    mtd[str(data.data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    # Dark current
    if dark_current is not None and dark_current.data is not None:
        subtract_dark_current(output_workspace, dark_current.data)

    # Normalization
    if flux_method is not None:
        kw = dict(method=flux_method, output_workspace=output_workspace)
        if flux_method == 'monitor':
            kw['monitor_workspace'] = data.monitor
        normalize_by_flux(output_workspace, flux, **kw)

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
                                 theta_deppendent_transmission=True,
                                 dark_current=None,
                                 flux_method=None,    # normalization (time/monitor/proton charge)
                                 flux=None,           # file for flux
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
                                 keep_processed_workspaces=True):
    r"""
    This function provides full data processing for a single experimental configuration,
    starting from workspaces (no data loading is happening inside this function)

    Parameters
    ----------
    sample_ws_raw: namedtuple
        (~mantid.dataobjects.Workspace2D, ~mantid.dataobjects.Workspace2D)
        raw data histogram workspace and monitor
    sample_trans_ws:  ~mantid.dataobjects.Workspace2D
        optional histogram workspace for sample transmission (already prepared)
    sample_trans_value: float
        optional value for sample transmission
    bkg_ws_raw: ~mantid.dataobjects.Workspace2D
        optional raw histogram workspace for background
    bkg_trans_ws: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for background transmission
    bkg_trans_value: float
        optional value for background transmission
    theta_deppendent_transmission: bool
        flag to apply angle dependent transmission
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
        flag to keep the processed background workspace

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    Sample workspace
    """
    if not output_workspace:
        output_workspace = output_suffix + '_sample'

    # create a common configuration for prepare data
    prepare_data_conf = {'dark_current': dark_current,
                         'flux_method': flux_method,
                         'flux': flux,
                         'mask_ws': mask_ws,
                         'mask_panel': mask_panel,
                         'mask_btp': mask_btp,
                         'solid_angle': solid_angle,
                         'sensitivity_workspace': sensitivity_workspace}

    # process sample
    sample_ws = prepare_data_workspaces(sample_ws_raw,
                                        output_workspace=output_workspace,
                                        **prepare_data_conf)
    # apply transmission to the sample
    if sample_trans_ws or sample_trans_value:
        sample_ws = apply_transmission_correction(sample_ws,
                                                  trans_workspace=sample_trans_ws,
                                                  trans_value=sample_trans_value,
                                                  theta_dependent=theta_deppendent_transmission,
                                                  output_workspace=output_workspace)

    # process background, if not already processed
    if bkg_ws_raw.data:
        bkgd_ws_name = output_suffix + '_background'
        if not registered_workspace(bkgd_ws_name):
            bkgd_ws = prepare_data_workspaces(bkg_ws_raw,
                                              output_workspace=bkgd_ws_name,
                                              **prepare_data_conf)
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

    # finalize with absolute scale and thickness
    sample_ws = normalize_by_thickness(sample_ws, thickness)

    # standard method assumes absolute scale from outside
    if absolute_scale_method == 'direct_beam':
        raise NotImplementedError("This feature is not yet implemented for EQSANS")
    else:
        sample_ws *= absolute_scale

    return mtd[output_workspace], sample_ws


def reduce_single_configuration(loaded_ws, reduction_input, prefix=''):
    flux_method = reduction_input["configuration"]["normalization"]
    if flux_method == "Monitor":
        flux_method = 'monitor'
        flux = reduction_input["configuration"]["fluxMonitorRatioFile"]
    elif flux_method == "Total charge":
        flux_method = 'proton charge'
        flux = reduction_input["configuration"]["beamFluxFileName"]
    elif flux_method == "Time":
        flux_method = 'time'
        flux = 'duration'
    else:
        flux_method = None
        flux = None

    solid_angle = reduction_input["configuration"]["useSolidAngleCorrection"]

    try:
        transmission_radius = float(reduction_input["configuration"]["mmRadiusForTransmission"])
    except ValueError:
        transmission_radius = None

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
    except (ValueError, KeyError):
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
    even_decades = reduction_input["configuration"]["LogQBinsEvenDecade"]
    nbins_main = int(reduction_input["configuration"]["numQBins"])
    outputFilename = reduction_input["outputFilename"]
    weighted_errors = reduction_input["configuration"]["useErrorWeighting"]
    try:
        qmin = float(reduction_input["configuration"]["Qmin"])
    except ValueError:
        qmin = None
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

    # automatically determine wedge binning if it wasn't explicitly set
    autoWedgeOpts = {}
    symmetric_wedges = True
    if bin1d_type == 'wedge':
        if wedges_min.size == 0:
            try:
                autoWedgeOpts = {'q_min': float(reduction_input['configuration']['autoWedgeQmin']),
                                 'q_delta': float(reduction_input['configuration']['autoWedgeQdelta']),
                                 'q_max': float(reduction_input['configuration']['autoWedgeQmax']),
                                 'azimuthal_delta': float(reduction_input['configuration']['autoWedgeAzimuthalDelta']),
                                 'peak_width': float(reduction_input['configuration'].get('autoWedgePeakWidth', 0.25)),
                                 'background_width': float(reduction_input['configuration']
                                                           .get('autoWedgeBackgroundWidth', 1.5)),
                                 'signal_to_noise_min': float(reduction_input['configuration']
                                                              .get('autoSignalToNoiseMin', 2.))}
                # auto-aniso returns all of the wedges
                symmetric_wedges = False
            except (KeyError, ValueError) as e:
                raise RuntimeError('Selected 1DQbinType="wedge", must either specify wedge angles or parameters '
                                   'for automatically determining them') from e

    # empty beam transmission workspace
    if loaded_ws.empty.data is not None:
        empty_trans_ws_name = f'{prefix}_empty'
        empty_trans_ws = prepare_data_workspaces(loaded_ws.empty,
                                                 flux_method=flux_method,
                                                 flux=flux,
                                                 solid_angle=False,
                                                 sensitivity_ws=loaded_ws.sensitivity,
                                                 output_workspace=empty_trans_ws_name)
    else:
        empty_trans_ws = None

    # background transmission
    background_transmission_dict = {}
    background_transmission_raw_dict = {}
    if loaded_ws.background_transmission.data is not None and empty_trans_ws is not None:
        bkgd_trans_ws_name = f'{prefix}_bkgd_trans'
        bkgd_trans_ws_processed = prepare_data_workspaces(loaded_ws.background_transmission,
                                                          flux_method=flux_method,
                                                          flux=flux,
                                                          solid_angle=False,
                                                          sensitivity_ws=loaded_ws.sensitivity,
                                                          output_workspace=bkgd_trans_ws_name)
        bkgd_trans_ws = calculate_transmission(bkgd_trans_ws_processed, empty_trans_ws,
                                               radius=transmission_radius, radius_unit="mm")
        print('Background transmission =', bkgd_trans_ws.extractY()[0, 0])

        bkgd_trans = reduction_input["background"]["transmission"]["runNumber"].strip()
        bk_tr_fn = os.path.join(output_dir, f'{outputFilename}_bkgd_{bkgd_trans}_trans.txt')
        SaveAscii(bkgd_trans_ws, Filename=bk_tr_fn)
        background_transmission_dict['value'] = bkgd_trans_ws.extractY()
        background_transmission_dict['error'] = bkgd_trans_ws.extractE()
        background_transmission_dict['wavelengths'] = bkgd_trans_ws.extractX()
        bkgd_trans_raw_ws = calculate_transmission(bkgd_trans_ws_processed, empty_trans_ws,
                                                   radius=transmission_radius, radius_unit="mm",
                                                   fit_function='')
        bk_tr_raw_fn = os.path.join(output_dir, f'{outputFilename}_bkgd_{bkgd_trans}_raw_trans.txt')
        SaveAscii(bkgd_trans_raw_ws, Filename=bk_tr_raw_fn)
        background_transmission_raw_dict['value'] = bkgd_trans_raw_ws.extractY()
        background_transmission_raw_dict['error'] = bkgd_trans_raw_ws.extractE()
        background_transmission_raw_dict['wavelengths'] = bkgd_trans_raw_ws.extractX()
    else:
        bkgd_trans_ws = None

    # sample transmission
    sample_transmission_dict = {}
    sample_transmission_raw_dict = {}
    if loaded_ws.sample_transmission.data is not None and empty_trans_ws is not None:
        sample_trans_ws_name = f'{prefix}_sample_trans'
        sample_trans_ws_processed = prepare_data_workspaces(loaded_ws.sample_transmission,
                                                            flux_method=flux_method,
                                                            flux=flux,
                                                            solid_angle=False,
                                                            sensitivity_ws=loaded_ws.sensitivity,
                                                            output_workspace=sample_trans_ws_name)
        sample_trans_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
                                                 radius=transmission_radius, radius_unit="mm")
        print('Sample transmission =', sample_trans_ws.extractY()[0, 0])

        tr_fn = os.path.join(output_dir, f'{outputFilename}_trans.txt')
        SaveAscii(sample_trans_ws, Filename=tr_fn)
        sample_transmission_dict['value'] = sample_trans_ws.extractY()
        sample_transmission_dict['error'] = sample_trans_ws.extractE()
        sample_transmission_dict['wavelengths'] = sample_trans_ws.extractX()

        sample_trans_raw_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
                                                     radius=transmission_radius, radius_unit="mm",
                                                     fit_function='')
        raw_tr_fn = os.path.join(output_dir, f'{outputFilename}_raw_trans.txt')
        SaveAscii(sample_trans_raw_ws, Filename=raw_tr_fn)
        sample_transmission_raw_dict['value'] = sample_trans_raw_ws.extractY()
        sample_transmission_raw_dict['error'] = sample_trans_raw_ws.extractE()
        sample_transmission_raw_dict['wavelengths'] = sample_trans_raw_ws.extractX()

    else:
        sample_trans_ws = None

    output = []
    detectordata = {}
    for i, raw_sample_ws in enumerate(loaded_ws.sample):
        name = "frame_{}".format(i+1)
        if len(loaded_ws.sample) > 1:
            output_suffix = f'_{i}'

        processed_data_main, sample_ws = process_single_configuration(raw_sample_ws,
                                                                      sample_trans_ws=sample_trans_ws,
                                                                      sample_trans_value=sample_trans_value,
                                                                      bkg_ws_raw=loaded_ws.background,
                                                                      bkg_trans_ws=bkgd_trans_ws,
                                                                      bkg_trans_value=bkg_trans_value,
                                                                      theta_deppendent_transmission=theta_deppendent_transmission,  # noqa E502
                                                                      dark_current=loaded_ws.dark_current,
                                                                      flux_method=flux_method,
                                                                      flux=flux,
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
        # Save nexus processed
        filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}.nxs')
        SaveNexus(processed_data_main, Filename=filename)
        # binning
        iq1d_main_in = convert_to_q(processed_data_main, mode='scalar')
        iq2d_main_in = convert_to_q(processed_data_main, mode='azimuthal')
        if bool(autoWedgeOpts):  # determine wedges automatically from the main detector
            wedges = getWedgeSelection(iq2d_main_in, **autoWedgeOpts)
        iq1d_main_in_fr = split_by_frame(processed_data_main, iq1d_main_in)
        iq2d_main_in_fr = split_by_frame(processed_data_main, iq2d_main_in)
        n_wl_frames = len(iq2d_main_in_fr)
        fr_label = ''
        for wl_frame in range(n_wl_frames):
            if n_wl_frames > 1:
                fr_label = f'_wl_frame_{wl_frame}'
            iq2d_main_out, iq1d_main_out = bin_all(iq2d_main_in_fr[wl_frame], iq1d_main_in_fr[wl_frame],
                                                   nxbins_main, nybins_main, nbins_main,
                                                   bin1d_type=bin1d_type, log_scale=log_binning,
                                                   even_decade=even_decades, qmin=qmin, qmax=qmax,
                                                   annular_angle_bin=annular_bin, wedges=wedges,
                                                   symmetric_wedges=symmetric_wedges,
                                                   error_weighted=weighted_errors)

            detectordata[name+fr_label] = {'iq': iq1d_main_out,
                                           'iqxqy': iq2d_main_out}

            # save ASCII files
            filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}{fr_label}_Iqxqy.dat')
            save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_main_out)

            for j in range(len(iq1d_main_out)):
                add_suffix = ""
                if len(iq1d_main_out) > 1:
                    add_suffix = f'_wedge_{j}'
                add_suffix += fr_label
                ascii_1D_filename = os.path.join(output_dir,
                                                 f'{outputFilename}{output_suffix}{add_suffix}_Iq.dat')
                save_iqmod(iq1d_main_out[j], ascii_1D_filename)

            IofQ_output = namedtuple('IofQ_output', ['I2D_main', 'I1D_main'])
            current_output = IofQ_output(I2D_main=iq2d_main_out,
                                         I1D_main=iq1d_main_out)
            output.append(current_output)

    # create reduction log
    filename = os.path.join(reduction_input["configuration"]["outputDir"], outputFilename + '_reduction_log.hdf')
    starttime = datetime.now().isoformat()
    pythonfile = __file__
    reductionparams = {'data': copy.deepcopy(reduction_input),
                       'filename': 'internal_file'}
    specialparameters = {'beam_center': {'x': 'not implemented yet',
                                         'y': 'not implemented yet',
                                         },
                         'sample_transmission': sample_transmission_dict,
                         'sample_transmission_raw': sample_transmission_raw_dict,
                         'background_transmission': background_transmission_dict,
                         'background_transmission_raw': background_transmission_raw_dict,
                         }

    samplelogs = {'main': SampleLogs(sample_ws)}
    drtsans.savereductionlog(filename=filename,
                             detectordata=detectordata,
                             reductionparams=reductionparams,
                             pythonfile=pythonfile,
                             starttime=starttime,
                             specialparameters=specialparameters,
                             samplelogs=samplelogs,
                             )

    return output


def plot_reduction_output(reduction_output, reduction_input, imshow_kwargs=None):
    output_dir = reduction_input["configuration"]["outputDir"]
    outputFilename = reduction_input["outputFilename"]
    output_suffix = ''

    bin1d_type = reduction_input["configuration"]["1DQbinType"]

    if imshow_kwargs is None:
        imshow_kwargs = {}
    for i, out in enumerate(reduction_output):
        if len(reduction_output) > 1:
            output_suffix = f'_{i}'

        if bin1d_type == 'wedge':
            wedges_min = np.fromstring(reduction_input["configuration"]["WedgeMinAngles"], sep=',')
            wedges_max = np.fromstring(reduction_input["configuration"]["WedgeMaxAngles"], sep=',')
            wedges = list(zip(wedges_min, wedges_max))
        else:
            wedges = None

        qmin = qmax = None
        if bin1d_type == 'wedge' or bin1d_type == 'annular':
            try:
                qmin = float(reduction_input["configuration"]["Qmin"])
            except ValueError:
                pass
            try:
                qmax = float(reduction_input["configuration"]["Qmax"])
            except ValueError:
                pass

        filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}_Iqxqy.png')
        plot_IQazimuthal(out.I2D_main, filename, backend='mpl',
                         imshow_kwargs=imshow_kwargs, title='Main',
                         wedges=wedges, qmin=qmin, qmax=qmax)

        for j in range(len(out.I1D_main)):
            add_suffix = ""
            if len(out.I1D_main) > 1:
                add_suffix = f'_wedge_{j}'
            filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}{add_suffix}_Iq.png')
            plot_IQmod([out.I1D_main[j]], filename, loglog=True,
                       backend='mpl', errorbar_kwargs={'label': 'main'})


def apply_solid_angle_correction(input_workspace):
    """Apply solid angle correction. This uses :func:`drtsans.solid_angle_correction`."""
    return solid_angle_correction(input_workspace,
                                  detector_type='VerticalTube')


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 center_x=None, center_y=None,
                 dark_current=None,
                 flux_method=None, flux=None,
                 mask=None, mask_panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None, sensitivity_workspace=None,
                 sample_aperture_diameter=None, sample_thickness=None,
                 source_aperture_diameter=None,
                 pixel_size_x=None, pixel_size_y=None,
                 output_workspace=None, output_suffix=''):
    r"""
    Load an EQSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    detector_offset: float
        Additional translation of the detector along Z-axis, in mili-meters.
    sample_offset: float
        Additional translation of the sample along the Z-axis, in mili-meters.
    bin_width: float
        Bin width for the output workspace, in Angstroms.
    low_tof_clip: float
        Ignore events with a time-of-flight (TOF) smaller than the minimal
        TOF plus this quantity.
    high_tof_clip: float
        Ignore events with a time-of-flight (TOF) bigger than the maximal
        TOF minus this quantity.
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'proton charge',
        'monitor', or 'time'.
    flux: str
        if ``flux_method`` is proton charge, then path to file containing the
        wavelength distribution of the neutron flux. If ``flux method`` is
        monitor, then path to file containing the flux-to-monitor ratios.
        if ``flux_method`` is time, then pass one log entry name such
        as ``duration`` or leave it as :py:obj:`None` for automatic log search.
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
    sample_aperture_diameter: float, None
        sample aperture diameter in unit mm
    sample_thickness: None, float
        sample thickness in unit cm
    source_aperture_diameter: float, None
        source aperture diameter in unit meter
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
    # First, load the event stream data into a workspace
    # The output_workspace name is for the Mantid workspace
    workspaces = load_events_and_histogram(data,
                                           detector_offset=detector_offset,
                                           sample_offset=sample_offset,
                                           output_workspace=output_workspace, output_suffix=output_suffix,
                                           center_x=center_x, center_y=center_y,
                                           bin_width=bin_width,
                                           low_tof_clip=low_tof_clip, high_tof_clip=high_tof_clip,
                                           keep_events=(dark_current is None),
                                           monitors=(flux_method == 'monitor'))

    output_workspace = workspaces.data

    # Next, we subtract dark current, if it exists.
    # Note that the function handles the normalization internally.
    if dark_current is not None:
        output_workspace = subtract_dark_current(output_workspace, dark_current)

    # The solid angle is corrected for next
    if solid_angle:
        if solid_angle is True:
            output_workspace = apply_solid_angle_correction(output_workspace)
        else:  # assume the solid_angle parameter is a workspace
            output_workspace = apply_solid_angle_correction(output_workspace, solid_angle_ws=solid_angle)

    # Interestingly, this is the only use of the btp dictionary.
    # The BTP stands for banks, tubes and pixels - it is a Mantid thing.
    apply_mask(output_workspace, panel=mask_panel, mask=mask, **btp)  # returns the mask

    # Correct for the detector sensitivity (the per pixel relative response)
    if sensitivity_file_path is not None or sensitivity_workspace is not None:
        kw = dict(sensitivity_filename=sensitivity_file_path, sensitivity_workspace=sensitivity_workspace)
        output_workspace = apply_sensitivity_correction(output_workspace, **kw)

    # We can perform the desired normalization here.
    if flux_method is not None:
        kw = dict(method=flux_method)
        if flux_method == 'monitor':
            kw['monitor_workspace'] = str(workspaces.monitor)
        output_workspace = normalize_by_flux(output_workspace, flux, **kw)

    # Overwrite meta data
    set_meta_data(output_workspace, None, None,
                  sample_offset,
                  sample_aperture_diameter, sample_thickness,
                  source_aperture_diameter,
                  pixel_size_x, pixel_size_y)

    if isinstance(output_workspace, str):
        return mtd[output_workspace]  # shouldn't happen
    else:
        return output_workspace
