""" BIOSANS API """
import copy
from datetime import datetime
import os
import numpy as np
import ast
from collections import namedtuple

from mantid.simpleapi import mtd, MaskDetectors

import drtsans
from drtsans import getWedgeSelection
from drtsans.path import registered_workspace
from drtsans.sensitivity import apply_sensitivity_correction, load_sensitivity_workspace
from drtsans.instruments import extract_run_number
from drtsans.samplelogs import SampleLogs
from drtsans.settings import namedtuplefy
from drtsans.plots import plot_IQmod, plot_IQazimuthal
from drtsans import subtract_background
from drtsans.reductionlog import savereductionlog
from drtsans.mono import biosans
from drtsans.mono.biosans import solid_angle_correction
from drtsans.mask_utils import apply_mask, load_mask
from drtsans.mono.load import load_events, transform_to_wavelength, set_init_uncertainties
from drtsans.mono.normalization import normalize_by_monitor, normalize_by_time
from drtsans.mono.dark_current import subtract_dark_current
from drtsans.mono.meta_data import set_meta_data
from drtsans.transmission import apply_transmission_correction, calculate_transmission
from drtsans.thickness_normalization import normalize_by_thickness
from drtsans.iq import bin_all
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D
from drtsans.path import allow_overwrite
from drtsans.dataobjects import IQmod
# from drtsans.mono.absolute_units import empty_beam_scaling
# from drtsans.mono.gpsans.attenuation import attenuation_factor

# Functions exposed to the general user (public) API
__all__ = ['prepare_data', 'load_all_files', 'plot_reduction_output',
           'prepare_data_workspaces', 'process_single_configuration',
           'reduce_single_configuration']


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
        blocked_beam = reduction_input["configuration"]["BlockBeamRunNumber"]
    else:
        blocked_beam = None

    # sample offsets, etc
    if load_params is None:
        load_params = {}
    try:
        wavelength = float(reduction_input["configuration"]["wavelength"])
        wavelength_spread_user = float(reduction_input["configuration"]["wavelengthSpread"])
    except ValueError:
        wavelength = wavelength_spread_user = None
    # if wavelength and wavelength_spread_user:
    #     # load_params["wavelength"] = wavelength
    #     # load_params["wavelengthSpread"] = wavelength_spread_user
    #     pass   # wavelength and wavelengthSpread not supported by LoadEventNexus at all
    # else:
    #     # reset user-overwriting wavelength and wavelength spread to None in case only 1 is specified
    #     wavelength = wavelength_spread_user = None

    if reduction_input["configuration"]["useDefaultMask"]:
        default_mask = [ast.literal_eval(mask_par) for mask_par in reduction_input["configuration"]["DefaultMask"]]
    else:
        default_mask = []

    path = f"/HFIR/{instrument_name}/IPTS-{ipts}/nexus"

    # check for time/log slicing
    timeslice = reduction_input["configuration"].get("timeslice")
    logslice = reduction_input["configuration"].get("logslice")

    if timeslice and logslice:
        raise ValueError("Can't do both time slicing and log slicing")

    if timeslice or logslice:
        if len(sample.split(',')) > 1:
            raise ValueError("Can't do slicing on summed data sets")

    # sample thickness
    try:
        # thickness is written to sample log if it is defined...
        # FIXME - thickness is used in reduce_configuration... - shall these 2 places more unified?
        thickness = float(reduction_input['thickness'])
    except ValueError:
        thickness = None

    # sample aperture diameter in mm
    try:
        sample_aperture_diameter = float(reduction_input['configuration']['sampleApertureSize'])
    except ValueError:
        sample_aperture_diameter = None
    except KeyError:
        raise KeyError('Please add "sampleApertureSize" under top-level section in JSON file')

    # source aperture diameter in mm
    try:
        source_aperture_diameter = float(reduction_input['configuration']['sourceApertureDiameter'])
    except ValueError:
        source_aperture_diameter = None
    except KeyError:
        raise KeyError('Please add "sourceApertureDiameter" under top-level section in JSON file')

    print(">>>>>>> reduction_input")
    print(reduction_input)

    # special loading case for sample to allow the slicing options
    logslice_data_dict = {}
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
            biosans.load_and_split(filename, output_workspace=ws_name,
                                   time_interval=timesliceinterval,
                                   log_name=logslicename, log_value_interval=logsliceinterval,
                                   **load_params)

            for _w in mtd[ws_name]:
                # Overwrite meta data
                set_meta_data(str(_w),
                              wave_length=wavelength,
                              wavelength_spread=wavelength_spread_user,
                              sample_thickness=thickness,
                              sample_aperture_diameter=sample_aperture_diameter,
                              source_aperture_diameter=source_aperture_diameter,
                              pixel_size_x=None,
                              pixel_size_y=None)
                # Transform X-axis to wave length with spread
                _w = transform_to_wavelength(_w)
                _w = set_init_uncertainties(_w)
                for btp_params in default_mask:
                    apply_mask(_w, **btp_params)

            for n in range(mtd[ws_name].getNumberOfEntries()):
                samplelogs = SampleLogs(mtd[ws_name].getItem(n))
                logslice_data_dict[str(n)] = {'data': list(samplelogs[logslicename].value),
                                              'units': samplelogs[logslicename].units,
                                              'name': logslicename}
    else:
        ws_name = f'{prefix}_{instrument_name}_{sample}_raw_histo'
        if not registered_workspace(ws_name):
            filename = ','.join(f"{path}/{instrument_name}_{run.strip()}.nxs.h5" for run in sample.split(','))
            print(f"Loading filename {filename}")
            biosans.load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
            # Overwrite meta data
            set_meta_data(ws_name,
                          wave_length=wavelength,
                          wavelength_spread=wavelength_spread_user,
                          sample_thickness=thickness,
                          sample_aperture_diameter=sample_aperture_diameter,
                          source_aperture_diameter=source_aperture_diameter,
                          pixel_size_x=None,
                          pixel_size_y=None)
            # Re-transform to wave length if overwriting values are specified
            if wavelength and wavelength_spread_user:
                transform_to_wavelength(ws_name)
            # Apply mask
            for btp_params in default_mask:
                apply_mask(ws_name, **btp_params)

    reduction_input["logslice_data"] = logslice_data_dict

    # load all other files
    for run_number in [center, bkgd, empty, sample_trans, bkgd_trans, blocked_beam]:
        if run_number:
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                filename = ','.join(f"{path}/{instrument_name}_{run.strip()}.nxs.h5" for run in run_number.split(','))
                print(f"Loading filename {filename}")
                biosans.load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)

    # do the same for dark current if exists
    dark_current_main = None
    dark_current_wing = None
    if reduction_input["configuration"]["useDarkFileName"]:
        dark_current_file_main = reduction_input["configuration"]["darkMainFileName"]
        dark_current_file_wing = reduction_input["configuration"]["darkWingFileName"]
        if dark_current_file_main and dark_current_file_wing:
            run_number = extract_run_number(dark_current_file_main)
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                print(f"Loading filename {dark_current_file_main}")
                dark_current_main = biosans.load_events_and_histogram(dark_current_file_main,
                                                                      output_workspace=ws_name,
                                                                      **load_params)
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)
            else:
                dark_current_main = mtd[ws_name]
            run_number = extract_run_number(dark_current_file_wing)
            ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
            if not registered_workspace(ws_name):
                print(f"Loading filename {dark_current_file_wing}")
                dark_current_wing = biosans.load_events_and_histogram(dark_current_file_wing,
                                                                      output_workspace=ws_name,
                                                                      **load_params)
                for btp_params in default_mask:
                    apply_mask(ws_name, **btp_params)
            else:
                dark_current_wing = mtd[ws_name]

    # load required processed_files
    sensitivity_main_ws_name = None
    sensitivity_wing_ws_name = None
    if reduction_input["configuration"]["useSensitivityFileName"]:
        flood_file_main = reduction_input["configuration"]["sensitivityMainFileName"]
        flood_file_wing = reduction_input["configuration"]["sensitivityWingFileName"]
        if flood_file_main and flood_file_wing:
            sensitivity_main_ws_name = f'{prefix}_main_sensitivity'
            sensitivity_wing_ws_name = f'{prefix}_wing_sensitivity'
            if not registered_workspace(sensitivity_main_ws_name):
                print(f"Loading filename {flood_file_main}")
                load_sensitivity_workspace(flood_file_main, output_workspace=sensitivity_main_ws_name)
            if not registered_workspace(sensitivity_wing_ws_name):
                print(f"Loading filename {flood_file_wing}")
                load_sensitivity_workspace(flood_file_wing, output_workspace=sensitivity_wing_ws_name)

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

    if registered_workspace(f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group'):
        raw_sample_ws = mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo_slice_group']
        raw_sample_ws_list = [w for w in raw_sample_ws]
    else:
        raw_sample_ws_list = [mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo']]
    raw_bkgd_ws = mtd[f'{prefix}_{instrument_name}_{bkgd}_raw_histo'] if bkgd else None
    raw_blocked_ws = mtd[f'{prefix}_{instrument_name}_{blocked_beam}_raw_histo'] if blocked_beam else None
    raw_center_ws = mtd[f'{prefix}_{instrument_name}_{center}_raw_histo']
    raw_empty_ws = mtd[f'{prefix}_{instrument_name}_{empty}_raw_histo'] if empty else None
    raw_sample_trans_ws = mtd[f'{prefix}_{instrument_name}_{sample_trans}_raw_histo'] if sample_trans else None
    raw_bkg_trans_ws = mtd[f'{prefix}_{instrument_name}_{bkgd_trans}_raw_histo'] if bkgd_trans else None
    sensitivity_main_ws = mtd[sensitivity_main_ws_name] if sensitivity_main_ws_name else None
    sensitivity_wing_ws = mtd[sensitivity_wing_ws_name] if sensitivity_wing_ws_name else None

    return dict(sample=raw_sample_ws_list,
                background=raw_bkgd_ws,
                center=raw_center_ws,
                empty=raw_empty_ws,
                sample_transmission=raw_sample_trans_ws,
                background_transmission=raw_bkg_trans_ws,
                blocked_beam=raw_blocked_ws,
                dark_current_main=dark_current_main,
                dark_current_wing=dark_current_wing,
                sensitivity_main=sensitivity_main_ws,
                sensitivity_wing=sensitivity_wing_ws,
                mask=mask_ws)


def prepare_data_workspaces(data,
                            center_x=None, center_y=None, center_y_wing=None,
                            dark_current=None,
                            flux_method=None,    # normalization (time/monitor)
                            mask_ws=None,        # apply a custom mask from workspace
                            mask_detector=None,  # main or wing
                            mask_panel=None,     # mask back or front panel
                            mask_btp=None,       # mask bank/tube/pixel
                            solid_angle=True,
                            sensitivity_workspace=None,
                            output_workspace=None):
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

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = str(data)
        output_workspace = output_workspace.replace('_raw_histo', '') + '_processed_histo'

    mtd[str(data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    if center_x is not None and center_y is not None and center_y_wing is not None:
        biosans.center_detector(output_workspace, center_x=center_x, center_y=center_y, center_y_wing=center_y_wing)

    # Dark current
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)

    # Normalization
    if str(flux_method).lower() == 'monitor':
        normalize_by_monitor(output_workspace)
    elif str(flux_method).lower() == 'time':
        normalize_by_time(output_workspace)

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
                                 center_y_wing=None,
                                 dark_current=None,
                                 flux_method=None,    # normalization (time/monitor)
                                 mask_ws=None,        # apply a custom mask from workspace
                                 mask_detector=None,
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
    center_y_wing: float
        y center for the wing
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
                         'center_y_wing': center_y_wing,
                         'dark_current': dark_current,
                         'flux_method': flux_method,
                         'mask_ws': mask_ws,
                         'mask_detector': mask_detector,
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
    transmission_dict = {}
    if sample_trans_ws or sample_trans_value:
        if sample_trans_value:
            transmission_dict = {'value': float(sample_trans_value),
                                 'error': ''}
        else:
            transmission_dict = {'value': sample_trans_ws.extractY(),
                                 'error': sample_trans_ws.extractE()}
        sample_ws = apply_transmission_correction(sample_ws,
                                                  trans_workspace=sample_trans_ws,
                                                  trans_value=sample_trans_value,
                                                  theta_dependent=theta_deppendent_transmission,
                                                  output_workspace=output_workspace)

    # process background, if not already processed
    background_transmission_dict = {}
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
                if bkg_trans_value:
                    background_transmission_dict = {'value': float(bkg_trans_value),
                                                    'error': ''}
                else:
                    background_transmission_dict = {'value': bkg_trans_ws.extractY(),
                                                    'error': bkg_trans_ws.extractE()}
                bkgd_ws = apply_transmission_correction(bkgd_ws,
                                                        trans_workspace=bkg_trans_ws,
                                                        trans_value=bkg_trans_value,
                                                        theta_dependent=theta_deppendent_transmission,
                                                        output_workspace=bkgd_ws_name)
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
    if absolute_scale_method == 'direct_beam':
        # try:
        #     empty = mtd[str(empty_beam_ws)]
        # except KeyError:
        #     raise ValueError(f"Could not find empty beam {str(empty_beam_ws)}")
        #
        # ac, ace = attenuation_factor(empty)
        # sample_ws = empty_beam_scaling(sample_ws,
        #                                empty,
        #                                beam_radius=beam_radius,
        #                                unit='mm',
        #                                attenuator_coefficient=ac,
        #                                attenuator_error=ace,
        #                                output_workspace=output_workspace)
        raise NotImplementedError('This method is not yet implemented for BIOSANS')
    else:
        sample_ws *= absolute_scale

    return mtd[str(sample_ws)], {'sample': transmission_dict,
                                 'background': background_transmission_dict}


def plot_reduction_output(reduction_output, reduction_input, loglog=True, imshow_kwargs=None):
    output_dir = reduction_input["configuration"]["outputDir"]
    outputFilename = reduction_input["outputFilename"]
    output_suffix = ''

    bin1d_type = reduction_input["configuration"]["1DQbinType"]

    if imshow_kwargs is None:
        imshow_kwargs = {}
    for i, out in enumerate(reduction_output):
        if len(reduction_output) > 1:
            output_suffix = f'_{i}'

        wedges = reduction_input["configuration"]["wedges"] if bin1d_type == 'wedge' else None
        symmetric_wedges = reduction_input["configuration"].get("symmetric_wedges", True)

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

        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D_main.png')
        plot_IQazimuthal(out.I2D_main, filename, backend='mpl',
                         imshow_kwargs=imshow_kwargs, title='Main',
                         wedges=wedges, symmetric_wedges=symmetric_wedges,
                         qmin=qmin, qmax=qmax)
        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D_wing.png')
        plot_IQazimuthal(out.I2D_wing, filename, backend='mpl',
                         imshow_kwargs=imshow_kwargs, title='Wing',
                         wedges=wedges, symmetric_wedges=symmetric_wedges,
                         qmin=qmin, qmax=qmax)
        for j in range(len(out.I1D_main)):
            add_suffix = ""
            if len(out.I1D_main) > 1:
                add_suffix = f'_wedge_{j}'
            filename = os.path.join(output_dir, '1D', f'{outputFilename}{output_suffix}_1D{add_suffix}.png')
            plot_IQmod([out.I1D_main[j], out.I1D_wing[j], out.I1D_combined[j]],
                       filename, loglog=loglog, backend='mpl', errorbar_kwargs={'label': 'main,wing,both'})

    # allow overwrite
    allow_overwrite(os.path.join(output_dir, '1D'))
    allow_overwrite(os.path.join(output_dir, '2D'))


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
        # thickness is written to sample log if it is defined...
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

    nxbins_main = int(reduction_input["configuration"]["numMainQxQyBins"])
    nybins_main = int(nxbins_main)
    nxbins_wing = int(reduction_input["configuration"]["numWingQxQyBins"])
    nybins_wing = int(nxbins_wing)
    bin1d_type = reduction_input["configuration"]["1DQbinType"]
    log_binning = reduction_input["configuration"]["QbinType"] == 'log'
    even_decades = reduction_input["configuration"].get("LogQBinsEvenDecade", False)
    decade_on_center = reduction_input["configuration"].get("LogQBinsDecadeCenter", False)
    try:
        nbins_main = int(reduction_input["configuration"].get("numMainQBins"))
    except (ValueError, KeyError, TypeError):
        nbins_main = None
    try:
        nbins_main_per_decade = int(reduction_input["configuration"].get("LogQBinsPerDecadeMain"))
    except (ValueError, KeyError, TypeError):
        nbins_main_per_decade = None
    try:
        nbins_wing = int(reduction_input["configuration"].get("numWingQBins"))
    except (ValueError, KeyError, TypeError):
        nbins_wing = None
    try:
        nbins_wing_per_decade = int(reduction_input["configuration"].get("LogQBinsPerDecadeWing"))
    except (ValueError, KeyError, TypeError):
        nbins_wing_per_decade = None
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

    xc, yc, yw = biosans.find_beam_center(loaded_ws.center)
    print("Center  =", xc, yc, yw)

    # empty beam transmission workspace
    if loaded_ws.empty is not None:
        empty_trans_ws_name = f'{prefix}_empty'
        empty_trans_ws = prepare_data_workspaces(loaded_ws.empty,
                                                 flux_method=flux_method,
                                                 mask_detector='wing_detector',
                                                 center_x=xc,
                                                 center_y=yc,
                                                 center_y_wing=yw,
                                                 solid_angle=False,
                                                 sensitivity_workspace=loaded_ws.sensitivity_main,
                                                 output_workspace=empty_trans_ws_name)
    else:
        empty_trans_ws = None

    # background transmission
    if loaded_ws.background_transmission is not None and empty_trans_ws is not None:
        bkgd_trans_ws_name = f'{prefix}_bkgd_trans'
        bkgd_trans_ws_processed = prepare_data_workspaces(loaded_ws.background_transmission,
                                                          flux_method=flux_method,
                                                          mask_detector='wing_detector',
                                                          center_x=xc,
                                                          center_y=yc,
                                                          center_y_wing=yw,
                                                          solid_angle=False,
                                                          sensitivity_workspace=loaded_ws.sensitivity_main,
                                                          output_workspace=bkgd_trans_ws_name)
        bkgd_trans_ws = calculate_transmission(bkgd_trans_ws_processed, empty_trans_ws,
                                               radius=transmission_radius, radius_unit="mm")
        print('Background transmission =', bkgd_trans_ws.extractY()[0, 0])
    else:
        bkgd_trans_ws = None

    # sample transmission
    if loaded_ws.sample_transmission is not None and empty_trans_ws is not None:
        sample_trans_ws_name = f'{prefix}_sample_trans'
        sample_trans_ws_processed = prepare_data_workspaces(loaded_ws.sample_transmission,
                                                            flux_method=flux_method,
                                                            mask_detector='wing_detector',
                                                            center_x=xc,
                                                            center_y=yc,
                                                            center_y_wing=yw,
                                                            solid_angle=False,
                                                            sensitivity_workspace=loaded_ws.sensitivity_main,
                                                            output_workspace=sample_trans_ws_name)
        sample_trans_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
                                                 radius=transmission_radius, radius_unit="mm")
        print('Sample transmission =', sample_trans_ws.extractY()[0, 0])
    else:
        sample_trans_ws = None

    output = []
    detectordata = {}
    for i, raw_sample_ws in enumerate(loaded_ws.sample):
        name = "_slice_{}".format(i+1)
        if len(loaded_ws.sample) > 1:
            output_suffix = f'_{i}'

        processed_data_main, trans_main = process_single_configuration(raw_sample_ws,
                                                                       sample_trans_ws=sample_trans_ws,
                                                                       sample_trans_value=sample_trans_value,
                                                                       bkg_ws_raw=loaded_ws.background,
                                                                       bkg_trans_ws=bkgd_trans_ws,
                                                                       bkg_trans_value=bkg_trans_value,
                                                                       blocked_ws_raw=loaded_ws.blocked_beam,
                                                                       theta_deppendent_transmission=theta_deppendent_transmission,  # noqa E502
                                                                       center_x=xc, center_y=yc, center_y_wing=yw,
                                                                       dark_current=loaded_ws.dark_current_main,
                                                                       flux_method=flux_method,
                                                                       mask_detector='wing_detector',
                                                                       mask_ws=loaded_ws.mask,
                                                                       mask_panel=mask_panel,
                                                                       solid_angle=solid_angle,
                                                                       sensitivity_workspace=loaded_ws.sensitivity_main,  # noqa E502
                                                                       output_workspace=f'processed_data_main_{i}',
                                                                       output_suffix=output_suffix,
                                                                       thickness=thickness,
                                                                       absolute_scale_method=absolute_scale_method,
                                                                       empty_beam_ws=empty_trans_ws,
                                                                       beam_radius=beam_radius,
                                                                       absolute_scale=absolute_scale,
                                                                       keep_processed_workspaces=False)
        processed_data_wing, trans_wing = process_single_configuration(raw_sample_ws,
                                                                       sample_trans_ws=sample_trans_ws,
                                                                       sample_trans_value=sample_trans_value,
                                                                       bkg_ws_raw=loaded_ws.background,
                                                                       bkg_trans_ws=bkgd_trans_ws,
                                                                       bkg_trans_value=bkg_trans_value,
                                                                       blocked_ws_raw=loaded_ws.blocked_beam,
                                                                       theta_deppendent_transmission=theta_deppendent_transmission,  # noqa E502
                                                                       center_x=xc, center_y=yc, center_y_wing=yw,
                                                                       dark_current=loaded_ws.dark_current_wing,
                                                                       flux_method=flux_method,
                                                                       mask_detector='detector1',
                                                                       mask_ws=loaded_ws.mask,
                                                                       mask_panel=mask_panel,
                                                                       solid_angle=solid_angle,
                                                                       sensitivity_workspace=loaded_ws.sensitivity_wing,  # noqa E502
                                                                       output_workspace=f'processed_data_wing_{i}',
                                                                       output_suffix=output_suffix,
                                                                       thickness=thickness,
                                                                       absolute_scale_method=absolute_scale_method,
                                                                       empty_beam_ws=empty_trans_ws,
                                                                       beam_radius=beam_radius,
                                                                       absolute_scale=absolute_scale,
                                                                       keep_processed_workspaces=False)

        # binning
        iq1d_main_in = biosans.convert_to_q(processed_data_main, mode='scalar')
        iq2d_main_in = biosans.convert_to_q(processed_data_main, mode='azimuthal')
        if bool(autoWedgeOpts):  # determine wedges automatically from the main detector
            wedges = getWedgeSelection(iq2d_main_in, **autoWedgeOpts)
            print('found wedge angles:')
            peak_wedge, back_wedge = wedges
            print('    peak:      ', peak_wedge)
            print('    background:', back_wedge)
            del peak_wedge, back_wedge

        # set the found wedge values to the reduction input, this will allow correct plotting
        reduction_input["configuration"]["wedges"] = wedges
        reduction_input["configuration"]["symmetric_wedges"] = symmetric_wedges

        iq2d_main_out, iq1d_main_out = bin_all(iq2d_main_in, iq1d_main_in, nxbins_main, nybins_main,
                                               n1dbins=nbins_main, n1dbins_per_decade=nbins_main_per_decade,
                                               decade_on_center=decade_on_center,
                                               bin1d_type=bin1d_type, log_scale=log_binning,
                                               even_decade=even_decades, qmin=qmin, qmax=qmax,
                                               annular_angle_bin=annular_bin, wedges=wedges,
                                               symmetric_wedges=symmetric_wedges,
                                               error_weighted=weighted_errors)
        iq1d_wing_in = biosans.convert_to_q(processed_data_wing, mode='scalar')
        iq2d_wing_in = biosans.convert_to_q(processed_data_wing, mode='azimuthal')
        iq2d_wing_out, iq1d_wing_out = bin_all(iq2d_wing_in, iq1d_wing_in, nxbins_wing, nybins_wing,
                                               n1dbins=nbins_wing, n1dbins_per_decade=nbins_wing_per_decade,
                                               decade_on_center=decade_on_center,
                                               bin1d_type=bin1d_type, log_scale=log_binning,
                                               even_decade=even_decades, qmin=qmin, qmax=qmax,
                                               annular_angle_bin=annular_bin, wedges=wedges,
                                               symmetric_wedges=symmetric_wedges,
                                               error_weighted=weighted_errors)

        # save ASCII files
        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D_main.dat')
        save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_main_out)
        filename = os.path.join(output_dir, '2D', f'{outputFilename}{output_suffix}_2D_wing.dat')
        save_ascii_binned_2D(filename, "I(Qx,Qy)", iq2d_wing_out)

        iq1d_combined_out = []
        for j in range(len(iq1d_main_out)):
            add_suffix = ""
            if len(iq1d_main_out) > 1:
                add_suffix = f'_wedge_{j}'
            ascii_1D_filename = os.path.join(output_dir, '1D',
                                             f'{outputFilename}{output_suffix}_1D_main{add_suffix}.txt')
            save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq1d_main_out[j])
            ascii_1D_filename = os.path.join(output_dir, '1D',
                                             f'{outputFilename}{output_suffix}_1D_wing{add_suffix}.txt')
            save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq1d_wing_out[j])
            try:
                OLT_Qmin = float(reduction_input["configuration"]["overlapStitchQmin"])
            except ValueError:
                OLT_Qmin = iq1d_wing_in.mod_q.min()
            try:
                OLT_Qmax = float(reduction_input["configuration"]["overlapStitchQmax"])
            except ValueError:
                OLT_Qmax = iq1d_main_in.mod_q.max()

            try:
                iq_output_both = biosans.stitch_profiles(profiles=[iq1d_main_out[j], iq1d_wing_out[j]],
                                                         overlaps=[OLT_Qmin, OLT_Qmax],
                                                         target_profile_index=0)

                ascii_1D_filename = os.path.join(output_dir, '1D',
                                                 f'{outputFilename}{output_suffix}_1D_both{add_suffix}.txt')
                save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq_output_both)
            except ZeroDivisionError:
                iq_output_both = IQmod(intensity=[], error=[], mod_q=[])
            iq1d_combined_out.append(iq_output_both)
        IofQ_output = namedtuple('IofQ_output', ['I2D_main', 'I2D_wing', 'I1D_main', 'I1D_wing', 'I1D_combined'])
        current_output = IofQ_output(I2D_main=iq2d_main_out,
                                     I2D_wing=iq2d_wing_out,
                                     I1D_main=iq1d_main_out,
                                     I1D_wing=iq1d_wing_out,
                                     I1D_combined=iq1d_combined_out)
        output.append(current_output)

        _inside_detectordata = {}
        if iq_output_both.intensity.size > 0:
            _inside_detectordata = {'combined': {'iq': [iq_output_both]}}
        else:
            _inside_detectordata = {}
        index = 0
        for _iq1d_main, _iq1d_wing, _iq2d_main, _iq2d_wing in zip(iq1d_main_out, iq1d_wing_out,
                                                                  [iq2d_main_out], [iq2d_wing_out]):
            _inside_detectordata["main_{}".format(index)] = {'iq': [_iq1d_main],
                                                             'iqxqy': _iq2d_main}
            _inside_detectordata["wing_{}".format(index)] = {'iq': [_iq1d_wing],
                                                             'iqxqy': _iq2d_wing}

        detectordata[name] = _inside_detectordata

    # save reduction log

    filename = os.path.join(reduction_input["configuration"]["outputDir"],
                            outputFilename + f'_reduction_log{output_suffix}.hdf')
    starttime = datetime.now().isoformat()
    # try:
    #     pythonfile = __file__
    # except NameError:
    #     pythonfile = "Launched from notebook"
    reductionparams = {'data': copy.deepcopy(reduction_input)}
    specialparameters = {'beam_center': {'x': xc,
                                         'y': yc,
                                         'y_wing': yw},
                         'sample_transmission': {'main': trans_main['sample'],
                                                 'wing': trans_wing['sample']},
                         'background_transmission': {'main': trans_main['background'],
                                                     'wing': trans_wing['background']}
                         }

    samplelogs = {'main': SampleLogs(processed_data_main),
                  'wing': SampleLogs(processed_data_wing)}
    logslice_data_dict = reduction_input['logslice_data']

    savereductionlog(filename=filename,
                     detectordata=detectordata,
                     reductionparams=reductionparams,
                     # pythonfile=pythonfile,
                     starttime=starttime,
                     specialparameters=specialparameters,
                     logslicedata=logslice_data_dict,
                     samplelogs=samplelogs,
                     )

    # change permissions to all files to allow overwrite
    allow_overwrite(reduction_input["configuration"]["outputDir"])
    allow_overwrite(os.path.join(reduction_input["configuration"]["outputDir"], '1D'))
    allow_overwrite(os.path.join(reduction_input["configuration"]["outputDir"], '2D'))

    return output


def prepare_data(data,
                 mask_detector=None,
                 detector_offset=0, sample_offset=0,
                 center_x=None, center_y=None, center_y_wing=None,
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
    Load a BIOSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    mask_detector: str
        Name of an instrument component to mask
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
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
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
    # TODO: missing detector_offset for wing detector
    ws = load_events(data, overwrite_instrument=True, output_workspace=output_workspace, output_suffix=output_suffix,
                     detector_offset=detector_offset, sample_offset=sample_offset)
    ws_name = str(ws)
    transform_to_wavelength(ws_name)
    set_init_uncertainties(ws_name)

    if center_x is not None and center_y is not None and center_y_wing is not None:
        biosans.center_detector(ws_name, center_x=center_x,
                                center_y=center_y,
                                center_y_wing=center_y_wing)

    # Mask either detector
    if mask_detector is not None:
        MaskDetectors(ws_name, ComponentList=mask_detector)

    # Dark current
    if dark_current is not None:
        if mtd.doesExist(str(dark_current)):
            dark_ws = mtd[str(dark_current)]
        else:
            dark_ws = load_events(dark_current, overwrite_instrument=True)
            dark_ws = transform_to_wavelength(dark_ws)
            dark_ws = set_init_uncertainties(dark_ws)
        subtract_dark_current(ws_name, dark_ws)

    # Normalization
    if str(flux_method).lower() == 'monitor':
        normalize_by_monitor(ws_name)
    elif str(flux_method).lower() == 'time':
        normalize_by_time(ws_name)

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
        drtsans.apply_sensitivity_correction(ws_name, sensitivity_filename=sensitivity_file_path,
                                             sensitivity_workspace=sensitivity_workspace)

    # Overwrite meta data
    set_meta_data(ws_name, wave_length, wavelength_spread,
                  sample_offset,
                  sample_aperture_diameter, sample_thickness,
                  source_aperture_diameter,
                  pixel_size_x, pixel_size_y)

    return mtd[ws_name]
