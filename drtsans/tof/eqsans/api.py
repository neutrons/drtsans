""" Top-level API for EQSANS """
from collections import namedtuple
import copy
from datetime import datetime
import os
import matplotlib.pyplot as plt
import numpy as np

from mantid.simpleapi import mtd, logger, SaveAscii, RebinToWorkspace, SaveNexus  # noqa E402
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
from drtsans.beam_finder import find_beam_center, fbc_options_json  # noqa E402
from drtsans.instruments import extract_run_number  # noqa E402
from drtsans.path import abspath, abspaths, registered_workspace  # noqa E402
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
from drtsans.path import allow_overwrite  # noqa E402
from drtsans.tof.eqsans.correction_api import CorrectionConfiguration
from drtsans.tof.eqsans.reduction_api import (prepare_data_workspaces, process_transmission)
from typing import Dict, Tuple


__all__ = ['apply_solid_angle_correction', 'subtract_background',
           'prepare_data', 'save_ascii_1D', 'save_xml_1D',
           'save_nist_dat', 'save_nexus', 'set_init_uncertainties',
           'load_all_files', 'prepare_data_workspaces',
           'pre_process_single_configuration', 'reduce_single_configuration',
           'plot_reduction_output']

IofQ_output = namedtuple('IofQ_output', ['I2D_main', 'I1D_main'])


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

    Workflow:
    1. parse reduction_input
    2. remove existing related workspaces with same run numbers
    3. process beam center
    -  output: load_params, reduction_input
    4. adjust pixel heights and widths
    -  output: load_params, reduction_input
    5. load and optionally slice sample runs
    6. load other runs: bkgd, empty, sample_trans, bkgd_trans

    Returned namedtuple:
        - sample, background, empty, sample_transmission, background_transmission: namedtuple(data[ws], monitor[ws])
        - dark_current, sensitivity, mask: workspace

    Returns
    -------
    namedtuple
        Named tuple including all loaded workspaces
    """
    reduction_config = reduction_input["configuration"]  # a handy shortcut to the configuration parameters dictionary

    instrument_name = reduction_input["instrumentName"]
    ipts = reduction_input["iptsNumber"]
    sample = reduction_input["sample"]["runNumber"]
    sample_trans = reduction_input["sample"]["transmission"]["runNumber"]
    bkgd = reduction_input["background"]["runNumber"]
    bkgd_trans = reduction_input["background"]["transmission"]["runNumber"]
    empty = reduction_input["emptyTransmission"]["runNumber"]
    center = reduction_input["beamCenter"]["runNumber"]

    # Remove existing workspaces, this is to guarantee that all the data is loaded correctly
    # In the future this should be made optional
    ws_to_remove = [f'{prefix}_{instrument_name}_{run_number}_raw_histo'
                    for run_number in (sample,
                                       bkgd,
                                       empty,
                                       sample_trans,
                                       bkgd_trans)]
    ws_to_remove.append(f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group')
    ws_to_remove.append(f'{prefix}_{instrument_name}_{center}_raw_events')
    ws_to_remove.append(f'{prefix}_sensitivity')
    ws_to_remove.append(f'{prefix}_mask')
    if reduction_config["darkFileName"]:
        run_number = extract_run_number(reduction_config["darkFileName"])
        ws_to_remove.append(f'{prefix}_{instrument_name}_{run_number}_raw_histo')
    for ws_name in ws_to_remove:
        if registered_workspace(ws_name):
            mtd.remove(ws_name)

    filenames = set()

    default_mask = None
    if reduction_config["useDefaultMask"]:
        configuration_file_parameters = _get_configuration_file_parameters(sample.split(',')[0].strip())
        default_mask = configuration_file_parameters['combined mask']

    # find the center first
    if center != "":
        # calculate beam center from center workspace
        center_ws_name = f'{prefix}_{instrument_name}_{center}_raw_events'
        if not registered_workspace(center_ws_name):
            center_filename = abspath(center, instrument=instrument_name, ipts=ipts)
            filenames.add(center_filename)
            load_events(center_filename,
                        pixel_calibration=reduction_config.get('usePixelCalibration', False),
                        output_workspace=center_ws_name)
            if reduction_config["useDefaultMask"]:
                apply_mask(center_ws_name, mask=default_mask)
        fbc_options = fbc_options_json(reduction_input)
        center_x, center_y, fit_results = find_beam_center(center_ws_name, **fbc_options)
        logger.notice(f"calculated center ({center_x}, {center_y})")
        beam_center_type = 'calculated'
    else:
        # use default EQSANS center
        # FIXME - it is better to have these hard code value defined out side of this method
        center_x = 0.025239
        center_y = 0.0170801
        logger.notice(f"use default center ({center_x}, {center_y})")
        beam_center_type = 'default'
    # set beam center
    reduction_input['beam_center'] = {'type': beam_center_type, 'x': center_x,
                                      'y': center_y, 'fit_results': fit_results}

    # update to 'load_params'
    if load_params is None:
        load_params = dict(center_x=center_x, center_y=center_y, keep_events=False)

    # Adjust pixel heights and widths
    load_params['pixel_calibration'] = reduction_config.get('usePixelCalibration', False)

    if reduction_config['detectorOffset'] is not None:
        load_params['detector_offset'] = reduction_config['detectorOffset']
    if reduction_config['sampleOffset'] is not None:
        load_params['sample_offset'] = reduction_config['sampleOffset']
    load_params['low_tof_clip'] = reduction_config["cutTOFmin"]
    load_params['high_tof_clip'] = reduction_config["cutTOFmax"]
    if reduction_config["wavelengthStep"] is not None:
        # account for wavelengthStepType
        step_type = 1
        if reduction_config["wavelengthStepType"] == "constant Delta lambda/lambda":
            step_type = -1
        load_params['bin_width'] = step_type * reduction_config["wavelengthStep"]
    load_params['monitors'] = reduction_config["normalization"] == "Monitor"

    # FIXME the issues with the monitor on EQSANS has not been fixed. Enable normalization by monitor (issue #538)
    if load_params['monitors']:
        raise RuntimeError('Normalization by monitor option will be enabled in a later drt-sans release')

    # check for time/log slicing
    timeslice,  logslice = reduction_config["useTimeSlice"], reduction_config["useLogSlice"]
    if timeslice or logslice:
        if len(sample.split(',')) > 1:
            raise ValueError("Can't do slicing on summed data sets")

    # Load (and optionally slice) sample runs
    # special loading case for sample to allow the slicing options
    logslice_data_dict = {}
    if timeslice or logslice:
        ws_name = f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group'
        if not registered_workspace(ws_name):
            filename = abspath(sample.strip(), instrument=instrument_name, ipts=ipts)
            print(f"Loading filename {filename}")
            if timeslice:
                timesliceinterval = reduction_config["timeSliceInterval"]
                logslicename = logsliceinterval = None
            elif logslice:
                timesliceinterval = None
                logslicename, logsliceinterval = reduction_config["logSliceName"], reduction_config["logSliceInterval"]
            filenames.add(filename)
            load_and_split(filename, output_workspace=ws_name,
                           time_interval=timesliceinterval,
                           log_name=logslicename, log_value_interval=logsliceinterval,
                           **load_params)
            for _w in mtd[ws_name]:
                if default_mask:
                    apply_mask(_w, mask=default_mask)

            if logslicename is not None:
                for n in range(mtd[ws_name].getNumberOfEntries()):
                    samplelogs = SampleLogs(mtd[ws_name].getItem(n))
                    logslice_data_dict[str(n)] = {'data': list(samplelogs[logslicename].value),
                                                  'units': samplelogs[logslicename].units,
                                                  'name': logslicename}
    else:
        ws_name = f'{prefix}_{instrument_name}_{sample}_raw_histo'
        if not registered_workspace(ws_name):
            filename = abspaths(sample.strip(), instrument=instrument_name, ipts=ipts)
            print(f"Loading filename {filename}")
            filenames.add(filename)
            load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
            if default_mask:
                apply_mask(ws_name, mask=default_mask)

    reduction_input["logslice_data"] = logslice_data_dict

    # load all other files without further processing
    # background, empty, sample transmission, background transmission
    for run_number in [bkgd, empty, sample_trans, bkgd_trans]:
        if run_number:
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                filename = abspaths(run_number.strip(), instrument=instrument_name, ipts=ipts)
                print(f"Loading filename {filename}")
                filenames.add(filename)
                load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
                if default_mask:
                    apply_mask(ws_name, mask=default_mask)

    # dark run (aka dark current run)
    dark_current_ws = None
    dark_current_mon_ws = None
    dark_current_file = reduction_config["darkFileName"]
    if dark_current_file is not None:
        run_number = extract_run_number(dark_current_file)
        ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
        if not registered_workspace(ws_name):
            dark_current_file = abspath(dark_current_file)
            print(f"Loading filename {dark_current_file}")
            filenames.add(dark_current_file)
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
    flood_file = reduction_input["configuration"]["sensitivityFileName"]
    if flood_file:
        sensitivity_ws_name = f'{prefix}_sensitivity'
        if not registered_workspace(sensitivity_ws_name):
            flood_file = abspath(flood_file)
            print(f"Loading filename {flood_file}")
            filenames.add(flood_file)
            load_sensitivity_workspace(flood_file, output_workspace=sensitivity_ws_name)
        sensitivity_ws = mtd[sensitivity_ws_name]

    mask_ws = None
    custom_mask_file = reduction_config["maskFileName"]
    if custom_mask_file is not None:
        mask_ws_name = f'{prefix}_mask'
        if not registered_workspace(mask_ws_name):
            custom_mask_file = abspath(custom_mask_file)
            print(f"Loading filename {custom_mask_file}")
            filenames.add(custom_mask_file)
            mask_ws = load_mask(custom_mask_file, output_workspace=mask_ws_name)
        else:
            mask_ws = mtd[mask_ws_name]

    # TODO load these files only once
    # beam_flux_ws = None
    # monitor_flux_ratio_ws = None

    sample_aperture_diameter = reduction_config["sampleApertureSize"]
    sample_thickness = reduction_input["sample"]["thickness"]
    smearing_pixel_size_x = reduction_config["smearingPixelSizeX"]
    smearing_pixel_size_y = reduction_config["smearingPixelSizeY"]

    # Sample workspace: set meta data
    for ws in sample_ws_list:
        set_meta_data(ws, wave_length=None, wavelength_spread=None,
                      sample_offset=load_params['sample_offset'],
                      sample_aperture_diameter=sample_aperture_diameter,
                      sample_thickness=sample_thickness,
                      source_aperture_diameter=None,
                      smearing_pixel_size_x=smearing_pixel_size_x,
                      smearing_pixel_size_y=smearing_pixel_size_y)

    print('FILE PATH, FILE SIZE:')
    total_size = 0
    for comma_separated_names in filenames:
        for name in comma_separated_names.split(','):
            try:
                file_size = os.path.getsize(name)
            except FileNotFoundError:
                hint = 'EQSANS_{}'.format(drtsans.instruments.extract_run_number(name))
                name = drtsans.path.abspath(hint, instrument='EQSANS')
                file_size = os.path.getsize(name)
            total_size += file_size
            print(name+',', '{:.2f} MiB'.format(file_size/1024**2))
    print('TOTAL: ', '{:.2f} MB'.format(total_size/1024**2))

    ws_mon_pair = namedtuple('ws_mon_pair', ['data', 'monitor'])

    loaded_ws_dict = dict(sample=[ws_mon_pair(data=ws, monitor=sample_mon_ws) for ws in sample_ws_list],
                          background=ws_mon_pair(data=background_ws, monitor=background_mon_ws),
                          empty=ws_mon_pair(data=empty_ws, monitor=empty_mon_ws),
                          sample_transmission=ws_mon_pair(data=sample_transmission_ws,
                                                          monitor=sample_transmission_mon_ws),
                          background_transmission=ws_mon_pair(data=background_transmission_ws,
                                                              monitor=background_transmission_mon_ws),
                          dark_current=ws_mon_pair(data=dark_current_ws, monitor=dark_current_mon_ws),
                          sensitivity=sensitivity_ws,
                          mask=mask_ws)

    return loaded_ws_dict


def pre_process_single_configuration(sample_ws_raw: namedtuple,
                                     sample_trans_ws=None,
                                     sample_trans_value=None,
                                     bkg_ws_raw=None,
                                     bkg_trans_ws=None,
                                     bkg_trans_value=None,
                                     theta_dependent_transmission=True,
                                     dark_current=None,
                                     flux_method=None,  # normalization (time/monitor/proton charge)
                                     flux=None,  # file for flux
                                     mask_ws=None,  # apply a custom mask from workspace
                                     mask_panel=None,  # mask back or front panel
                                     mask_btp=None,  # mask bank/tube/pixel
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
    theta_dependent_transmission: bool
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
        if sample_trans_ws:
            RebinToWorkspace(WorkspaceToRebin=sample_trans_ws,
                             WorkspaceToMatch=sample_ws,
                             OutputWorkspace=sample_trans_ws)
        sample_ws = apply_transmission_correction(sample_ws,
                                                  trans_workspace=sample_trans_ws,
                                                  trans_value=sample_trans_value,
                                                  theta_dependent=theta_dependent_transmission,
                                                  output_workspace=output_workspace)

    # process background, if not already processed
    if bkg_ws_raw.data:
        # process background run
        bkgd_ws_name = output_suffix + '_background'
        if not registered_workspace(bkgd_ws_name):
            bkgd_ws = prepare_data_workspaces(bkg_ws_raw,
                                              output_workspace=bkgd_ws_name,
                                              **prepare_data_conf)
            # apply transmission to background
            if bkg_trans_ws or bkg_trans_value:
                if bkg_trans_ws:
                    RebinToWorkspace(WorkspaceToRebin=bkg_trans_ws,
                                     WorkspaceToMatch=bkgd_ws,
                                     OutputWorkspace=bkg_trans_ws)
                bkgd_ws = apply_transmission_correction(bkgd_ws,
                                                        trans_workspace=bkg_trans_ws,
                                                        trans_value=bkg_trans_value,
                                                        theta_dependent=theta_dependent_transmission,
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

    return mtd[output_workspace]


# TODO 777 - better documentation of this function
# FIXME TODO Task 786 Remove temp_debug_weighting
# FIXME 777 Before closing: incoherence_correction_setup shall be set up inside reduce_single_configuration()
def reduce_single_configuration(loaded_ws: namedtuple,
                                reduction_input,
                                prefix='',
                                skip_nan=True,
                                incoherence_correction_setup=None):
    """Reduce samples from raw workspaces including
    1. prepare data
    1.

    This is the main entry point of reduction

    Input loaded workspaces as namedtuple:
      - sample, background, empty, sample_transmission, background_transmission: namedtuple(data[ws], monitor[ws])
      - dark_current, sensitivity, mask: workspace

    Parameters
    ----------
    loaded_ws: namedtuple
        loaded workspaces
    reduction_input: dict
        reduction configuration
    prefix
    skip_nan
    incoherence_correction_setup: CorrectionConfiguration, None
        incoherence/inelastic scattering correction configuration

    Returns
    -------
    ~list
        list of IofQ_output: ['I2D_main', 'I1D_main']

    """
    # Process reduction input: configuration and etc.
    reduction_config = reduction_input["configuration"]

    # Process inelastic/incoherent scattering correction configuration if user does not specify
    if incoherence_correction_setup is None:
        # backward compatibility
        # TODO 777 Task 1 - implement parse_correction_config() and unify with sections of codes in tests
        # incoherence_correction_setup = parse_correction_config(reduction_input)
        incoherence_correction_setup = CorrectionConfiguration(do_correction=False)

    # process: flux, monitor, proton charge, ...
    flux_method_translator = {'Monitor': 'monitor', 'Total charge': 'proton charge', 'Time': 'time'}
    flux_method = flux_method_translator.get(reduction_config["normalization"], None)

    flux_translator = {'Monitor': reduction_config["fluxMonitorRatioFile"],
                       'Total charge': reduction_config["beamFluxFileName"],
                       'Time': 'duration'}
    flux = flux_translator.get(reduction_config["normalization"], None)

    solid_angle = reduction_config["useSolidAngleCorrection"]
    transmission_radius = reduction_config["mmRadiusForTransmission"]
    sample_trans_value = reduction_input["sample"]["transmission"]["value"]
    bkg_trans_value = reduction_input["background"]["transmission"]["value"]
    theta_dependent_transmission = reduction_config["useThetaDepTransCorrection"]
    mask_panel = 'back' if reduction_config["useMaskBackTubes"] is True else None
    output_suffix = ''

    thickness = reduction_input["sample"]["thickness"]
    absolute_scale_method = reduction_config["absoluteScaleMethod"]
    beam_radius = None  # EQSANS doesn't use keyword DBScalingBeamRadius
    absolute_scale = reduction_config["StandardAbsoluteScale"]
    output_dir = reduction_config["outputDir"]
    # process binning
    # Note: option {even_decades = reduction_config["useLogQBinsEvenDecade"]} is removed
    nybins_main = nxbins_main = reduction_config["numQxQyBins"]
    bin1d_type = reduction_config["1DQbinType"]
    log_binning = reduction_config["QbinType"] == 'log'  # FIXME - note: fixed to log binning
    decade_on_center = reduction_config["useLogQBinsDecadeCenter"]
    nbins_main = reduction_config["numQBins"]
    nbins_main_per_decade = reduction_config["LogQBinsPerDecade"]
    outputFilename = reduction_input["outputFileName"]
    weighted_errors = reduction_config["useErrorWeighting"]
    qmin = reduction_config["Qmin"]
    qmax = reduction_config["Qmax"]
    annular_bin = reduction_config["AnnularAngleBin"]
    wedges_min = reduction_config["WedgeMinAngles"]
    wedges_max = reduction_config["WedgeMaxAngles"]
    wedges = None if wedges_min is None or wedges_max is None else list(zip(wedges_min, wedges_max))
    # set the found wedge values to the reduction input, this will allow correct plotting
    reduction_config["wedges"] = wedges
    reduction_config["symmetric_wedges"] = True

    # automatically determine wedge binning if it wasn't explicitly set
    autoWedgeOpts, symmetric_wedges = parse_auto_wedge_setup(reduction_config, bin1d_type, wedges_min)

    # Prepare empty beam transmission workspace
    if loaded_ws.empty.data is not None:
        empty_trans_ws_name = f'{prefix}_empty'
        empty_trans_ws = prepare_data_workspaces(loaded_ws.empty,
                                                 flux_method=flux_method,
                                                 flux=flux,
                                                 solid_angle=False,
                                                 sensitivity_workspace=loaded_ws.sensitivity,
                                                 output_workspace=empty_trans_ws_name)
    else:
        empty_trans_ws = None

    # Background transmission
    # specific output filename (base) for background trans
    bkgd_trans = reduction_input["background"]["transmission"]["runNumber"].strip()
    base_out_name = f'{outputFilename}_bkgd_{bkgd_trans}'

    # process transmission
    bkgd_returned = process_transmission(loaded_ws.background_transmission,
                                         empty_trans_ws,
                                         transmission_radius,
                                         loaded_ws.sensitivity,
                                         flux_method, flux, prefix, 'bkgd',
                                         output_dir, base_out_name)
    bkgd_trans_ws, background_transmission_dict, background_transmission_raw_dict = bkgd_returned

    # sample transmission
    sample_returned = process_transmission(loaded_ws.sample_transmission,
                                           empty_trans_ws,
                                           transmission_radius,
                                           loaded_ws.sensitivity,
                                           flux_method, flux, prefix, 'sample',
                                           output_dir, outputFilename)

    sample_trans_ws, sample_transmission_dict, sample_transmission_raw_dict = sample_returned

    # Form binning parameters
    # binning_par_dc = {'nxbins_main': nxbins_main,
    #                   'nybins_main': nybins_main,
    #                   'n1dbins': nbins_main,
    #                   'n1dbins_per_decade': nbins_main_per_decade,
    #                   'decade_on_center': decade_on_center,
    #                   'bin1d_type': bin1d_type,
    #                   'log_scale': log_binning,
    #                   'qmin': qmin,
    #                   'qmax': qmax,
    #                   'qxrange': None,
    #                   'qyrange': None}

    # binning_params = namedtuple('binning_setup', binning_par_dc)(**binning_par_dc)
    # binning_params = BinningSetup(**binning_par_dc)
    #

    # Define output data structure
    output = []
    detectordata = {}
    processed_data_main = None
    for i, raw_sample_ws in enumerate(loaded_ws.sample):
        name = "slice_{}".format(i+1)
        if len(loaded_ws.sample) > 1:
            output_suffix = f'_{i}'
        raw_name = f'EQSANS_{raw_sample_ws.data.getRunNumber()}'

        # process data without correction
        processed_data_main = pre_process_single_configuration(raw_sample_ws,
                                                               sample_trans_ws=sample_trans_ws,
                                                               sample_trans_value=sample_trans_value,
                                                               bkg_ws_raw=loaded_ws.background,
                                                               bkg_trans_ws=bkgd_trans_ws,
                                                               bkg_trans_value=bkg_trans_value,
                                                               theta_dependent_transmission=theta_dependent_transmission,  # noqa E502
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

        # Convert to Q
        # set up subpixel binning options  FIXME - it does not seem to work
        subpixel_kwargs = dict()
        if reduction_config['useSubpixels']:
            subpixel_kwargs = {'n_horizontal': reduction_config['subpixelsX'],
                               'n_vertical': reduction_config['subpixelsY']}
        # convert to Q
        iq1d_main_in = convert_to_q(processed_data_main, mode='scalar', **subpixel_kwargs)
        iq2d_main_in = convert_to_q(processed_data_main, mode='azimuthal', **subpixel_kwargs)

        # split to frames
        iq1d_main_in_fr = split_by_frame(processed_data_main, iq1d_main_in, verbose=True)
        iq2d_main_in_fr = split_by_frame(processed_data_main, iq2d_main_in, verbose=True)

        # Save nexus processed
        filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}.nxs')
        SaveNexus(processed_data_main, Filename=filename)
        print(f'SaveNexus to {filename}')

        # Work with wedges
        if bool(autoWedgeOpts):  # determine wedges automatically from the main detectora
            wedges = process_auto_wedge(autoWedgeOpts, iq2d_main_in, output_dir, reduction_config, symmetric_wedges)

        n_wl_frames = len(iq2d_main_in_fr)
        _inside_detectordata = {}

        # Process each frame separately
        user_qmin = qmin
        user_qmax = qmax

        for wl_frame in range(n_wl_frames):
            if n_wl_frames > 1:
                fr_log_label = f'_frame_{wl_frame}'
                fr_label = fr_log_label
            else:
                fr_log_label = f'frame'
                fr_label = ""

            assert iq1d_main_in_fr[wl_frame] is not None, 'Input I(Q)      main input cannot be None.'
            assert iq2d_main_in_fr[wl_frame] is not None, 'Input I(qx, qy) main input cannot be None.'

            # print(f'[DEBUG 777] Frame {wl_frame} '
            #       f'Before Q1D NaN = {np.where(np.isnan(iq1d_main_in_fr[wl_frame].intensity))[0].shape}; '
            #       f'Total data points = {iq1d_main_in_fr[wl_frame].intensity.shape}')
            # print(f'DEBUG 777] Is log binning: {log_binning}')
            # if True:
            #     qmin_fw = iq1d_main_in_fr[wl_frame].mod_q.min()
            #     qmax_fw = iq1d_main_in_fr[wl_frame].mod_q.max()
            #     print(f'[DEBUG 777] Frame {wl_frame} Q range: {qmin_fw}, {qmax_fw}')

            # Note 777: 2 step binning shall generate the same result as 1 step binning
            # FIXME 786 - Remove temp_debug_weighting after conceptual proved
            if incoherence_correction_setup.do_correction:  # or temp_debug_weighting:
                assert weighted_errors, 'Must using weighted error'

                print(f'[DEBUG 777-1A] Frame {wl_frame} '
                      f'Before Q1D NaN = {np.where(np.isnan(iq1d_main_in_fr[wl_frame].intensity))[0].shape}')

                # Define qmin and qmax for this frame
                if user_qmin is None:
                    qmin = iq1d_main_in_fr[wl_frame].mod_q.min()
                else:
                    qmin = user_qmin
                if user_qmax is None:
                    qmax = iq1d_main_in_fr[wl_frame].mod_q.max()
                else:
                    qmax = user_qmax
                print(f'[DEBUG 777-1] Frame {wl_frame} Q range: {qmin}, {qmax}')

                # Determine qxrange and qyrange for this frame
                qx_min = np.min(iq2d_main_in_fr[wl_frame].qx)
                qx_max = np.max(iq2d_main_in_fr[wl_frame].qx)
                qxrange = qx_min, qx_max
                print(f'[DEBUG 777-1] Frame {wl_frame} Qx range: {qx_min}, {qx_max}; User range: {qxrange}')

                qy_min = np.min(iq2d_main_in_fr[wl_frame].qy)
                qy_max = np.max(iq2d_main_in_fr[wl_frame].qy)
                qyrange = qy_min, qy_max
                print(f'[DEBUG 777-1] Frame {wl_frame} Qy range: {qy_min}, {qy_max}; User range: {qyrange}')

                # Bin I(Q1D, wl) and I(Q2D, wl) in Q and (Qx, Qy) space respectively but not wavelength
                iq2d_main_wl, iq1d_main_wl = bin_all(iq2d_main_in_fr[wl_frame], iq1d_main_in_fr[wl_frame],
                                                     nxbins_main, nybins_main, n1dbins=nbins_main,
                                                     n1dbins_per_decade=nbins_main_per_decade,
                                                     decade_on_center=decade_on_center,
                                                     bin1d_type=bin1d_type, log_scale=log_binning,
                                                     qmin=qmin, qmax=qmax,
                                                     qxrange=qxrange,
                                                     qyrange=qyrange,
                                                     annular_angle_bin=annular_bin, wedges=wedges,
                                                     symmetric_wedges=symmetric_wedges,
                                                     error_weighted=weighted_errors,
                                                     n_wavelength_bin=None)
                assert isinstance(iq1d_main_wl, list), f'Output I(Q) must be a list but not a {type(iq1d_main_wl)}'
                print(f'[DEBUG 777-1B] Middle Q1D NaN = {np.where(np.isnan(iq1d_main_wl[0].intensity))[0].shape}')

                if len(iq1d_main_wl) != 1:
                    raise NotImplementedError(f'Not expected that there are more than 1 IQmod main but '
                                              f'{len(iq1d_main_wl)}')

                # TODO 786 - Implement 1D
                # TODO 786 - Implement 1D correction
                from drtsans.tof.eqsans.correction_api import (do_inelastic_incoherence_correction_q1d,
                                                               do_inelastic_incoherence_correction_q2d)
                b_file_prefix = f'{raw_name}_frame_{wl_frame}'
                corrected_iq1d = do_inelastic_incoherence_correction_q1d(iq1d_main_wl[0],
                                                                         incoherence_correction_setup,
                                                                         b_file_prefix,
                                                                         output_dir)

                # TODO 787 - Implement 2D correction
                # TODO 787 - Implement 2D
                corrected_iq2d = do_inelastic_incoherence_correction_q2d(iq2d_main_wl,
                                                                         incoherence_correction_setup,
                                                                         b_file_prefix,
                                                                         output_dir)

                # Be finite
                finite_iq1d = corrected_iq1d.be_finite()
                finite_iq2d = corrected_iq2d.be_finite()
                # Bin binned I(Q1D, wl) and and binned I(Q2D, wl) in wavelength space
                assert len(iq1d_main_wl) == 1, f'It is assumed that output I(Q) list contains 1 I(Q)' \
                                               f' but not {len(iq1d_main_wl)}'
                iq2d_main_out, iq1d_main_out = bin_all(finite_iq2d, finite_iq1d,
                                                       nxbins_main, nybins_main, n1dbins=nbins_main,
                                                       n1dbins_per_decade=nbins_main_per_decade,
                                                       decade_on_center=decade_on_center,
                                                       bin1d_type=bin1d_type, log_scale=log_binning,
                                                       qmin=qmin, qmax=qmax,
                                                       qxrange=None,
                                                       qyrange=None,
                                                       annular_angle_bin=annular_bin, wedges=wedges,
                                                       symmetric_wedges=symmetric_wedges,
                                                       error_weighted=weighted_errors)

            else:
                # Not incoherence correction
                print(f'[DEBUG 777-2] Frame {wl_frame} '
                      f'Before Q1D NaN = {np.isnan(iq1d_main_in_fr[wl_frame].intensity).shape}')

                iq2d_main_out, iq1d_main_out = bin_all(iq2d_main_in_fr[wl_frame], iq1d_main_in_fr[wl_frame],
                                                       nxbins_main, nybins_main, n1dbins=nbins_main,
                                                       n1dbins_per_decade=nbins_main_per_decade,
                                                       decade_on_center=decade_on_center,
                                                       bin1d_type=bin1d_type, log_scale=log_binning,
                                                       qmin=qmin, qmax=qmax,
                                                       qxrange=None,
                                                       qyrange=None,
                                                       annular_angle_bin=annular_bin, wedges=wedges,
                                                       symmetric_wedges=symmetric_wedges,
                                                       error_weighted=weighted_errors)

            _inside_detectordata[fr_log_label] = {'iq': iq1d_main_out, 'iqxqy': iq2d_main_out}

            # save ASCII files
            filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}{fr_label}_Iqxqy.dat')
            # [#689] TODO FIXME - make save_ascii_binned_2D() back afterwards: iq2d_main_out cannot be None
            if iq2d_main_out:
                save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_main_out)

            for j in range(len(iq1d_main_out)):
                add_suffix = ""
                if len(iq1d_main_out) > 1:
                    add_suffix = f'_wedge_{j}'
                add_suffix += fr_label
                ascii_1D_filename = os.path.join(output_dir,
                                                 f'{outputFilename}{output_suffix}{add_suffix}_Iq.dat')
                save_iqmod(iq1d_main_out[j], ascii_1D_filename, skip_nan=skip_nan)

            current_output = IofQ_output(I2D_main=iq2d_main_out,
                                         I1D_main=iq1d_main_out)
            output.append(current_output)
        # END binning loop over frame

        detectordata[name] = _inside_detectordata
    # END reduction loop over sample workspaces

    # Save reduction log
    save_reduction_log(reduction_input, outputFilename, processed_data_main,
                       sample_transmission_dict, sample_transmission_raw_dict,
                       background_transmission_dict, background_transmission_raw_dict,
                       detectordata, output_dir)

    return output


def save_reduction_log(reduction_input, outputFilename, processed_data_main,
                       sample_transmission_dict, sample_transmission_raw_dict,
                       background_transmission_dict, background_transmission_raw_dict,
                       detectordata, output_dir):
    # TODO 784 Doc!

    # create reduction log
    filename = os.path.join(reduction_input["configuration"]["outputDir"],
                            outputFilename + f'_reduction_log.hdf')
    starttime = datetime.now().isoformat()
    # try:
    #     pythonfile = __file__
    # except NameError:
    #     pythonfile = "Launched from notebook"
    reductionparams = {'data': copy.deepcopy(reduction_input)}
    beam_center_dict = reduction_input['beam_center']
    specialparameters = {'beam_center': {'x': beam_center_dict['x'],
                                         'y': beam_center_dict['y'],
                                         'type': beam_center_dict['type']},
                         'fit_results': beam_center_dict['fit_results'],
                         'sample_transmission': sample_transmission_dict,
                         'sample_transmission_raw': sample_transmission_raw_dict,
                         'background_transmission': background_transmission_dict,
                         'background_transmission_raw': background_transmission_raw_dict}

    # [#689] TODO FIXME - Reincarnate this section!
    # FIXME - check original code.  processed_data_main outside a loop is a BUG!
    #  The correction workflow does not output processed data workspace yet!
    assert processed_data_main is not None
    samplelogs = {'main': SampleLogs(processed_data_main)}
    logslice_data_dict = reduction_input["logslice_data"]

    drtsans.savereductionlog(filename=filename,
                             detectordata=detectordata,
                             reductionparams=reductionparams,
                             starttime=starttime,
                             specialparameters=specialparameters,
                             logslicedata=logslice_data_dict,
                             samplelogs=samplelogs)

    # change permissions to all files to allow overwrite
    allow_overwrite(output_dir)


def process_auto_wedge(auto_wedge_setup: Dict, iq2d_input, output_dir: str, reduction_config: Dict,
                       symmetric_wedges):
    # TODO 784 Doc!
    logger.notice(f'Auto wedge options: {auto_wedge_setup}')
    auto_wedge_setup['debug_dir'] = output_dir
    wedges = getWedgeSelection(iq2d_input, **auto_wedge_setup)
    logger.notice(f'found wedge angles:\n'
                  f'              peak: {wedges[0]}\n'
                  f'        background: {wedges[1]}')
    # sanity check
    assert len(wedges) == 2, f'Auto-wedges {wedges} shall have 2 2-tuples'
    # set automated wedge to reduction configuration for correct plotting.
    # reduction_config is an in/out function argument
    reduction_config["wedges"] = wedges
    reduction_config["symmetric_wedges"] = symmetric_wedges

    return wedges


def parse_auto_wedge_setup(reduction_config: Dict, bin1d_type: str, wedges_min) -> Tuple[Dict, bool]:
    # TODO 784 Doc!
    autoWedgeOpts = {}
    symmetric_wedges = True
    if bin1d_type == 'wedge' and len(wedges_min) == 0:
        # the JSON validator "wedgesources" guarantees that the parameters to be collected are all non-empty
        autoWedgeOpts = {'q_min': reduction_config['autoWedgeQmin'],
                         'q_delta': reduction_config['autoWedgeQdelta'],
                         'q_max': reduction_config['autoWedgeQmax'],
                         'azimuthal_delta': reduction_config['autoWedgeAzimuthalDelta'],
                         'peak_width': reduction_config['autoWedgePeakWidth'],
                         'background_width': reduction_config['autoWedgeBackgroundWidth'],
                         'signal_to_noise_min': reduction_config['autoWedgeSignalToNoiseMin']}
        # auto-aniso returns all of the wedges
        symmetric_wedges = False

    return autoWedgeOpts, symmetric_wedges


def plot_reduction_output(reduction_output, reduction_input, imshow_kwargs=None):
    reduction_config = reduction_input['configuration']
    output_dir = reduction_config["outputDir"]
    outputFilename = reduction_input["outputFileName"]
    output_suffix = ''

    bin1d_type = reduction_config["1DQbinType"]

    if imshow_kwargs is None:
        imshow_kwargs = {}
    for i, out in enumerate(reduction_output):
        if len(reduction_output) > 1:
            output_suffix = f'_{i}'

        wedges = reduction_config["wedges"] if bin1d_type == 'wedge' else None
        symmetric_wedges = reduction_config.get("symmetric_wedges", True)

        qmin = reduction_config["Qmin"]
        qmax = reduction_config["Qmax"]

        filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}_Iqxqy.png')
        plot_IQazimuthal(out.I2D_main, filename, backend='mpl',
                         imshow_kwargs=imshow_kwargs, title='Main',
                         wedges=wedges, symmetric_wedges=symmetric_wedges,
                         qmin=qmin, qmax=qmax)
        plt.clf()
        for j in range(len(out.I1D_main)):
            add_suffix = ""
            if len(out.I1D_main) > 1:
                add_suffix = f'_wedge_{j}'
            filename = os.path.join(output_dir, f'{outputFilename}{output_suffix}{add_suffix}_Iq.png')
            plot_IQmod([out.I1D_main[j]], filename, loglog=True,
                       backend='mpl', errorbar_kwargs={'label': 'main'})
            plt.clf()
    plt.close()
    # change permissions to all files to allow overwrite
    allow_overwrite(output_dir)


def apply_solid_angle_correction(input_workspace):
    """Apply solid angle correction. This uses :func:`drtsans.solid_angle_correction`."""
    return solid_angle_correction(input_workspace,
                                  detector_type='VerticalTube')


def prepare_data(data,
                 pixel_calibration=False,
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
                 smearing_pixel_size_x=None, smearing_pixel_size_y=None,
                 output_workspace=None, output_suffix=''):
    r"""
    Load an EQSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    pixel_calibration: bool
        Adjust pixel heights and widths according to barscan and tube-width calibrations.
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
    smearing_pixel_size_x: float, None
        pixel size in x direction in unit as meter, only for Q-resolution calculation
    smearing_pixel_size_y: float, None
        pixel size in Y direction in unit as meter, only for Q-resolutio calculation

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
                                           pixel_calibration=pixel_calibration,
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
                  smearing_pixel_size_x, smearing_pixel_size_y)

    if isinstance(output_workspace, str):
        return mtd[output_workspace]  # shouldn't happen
    else:
        return output_workspace
