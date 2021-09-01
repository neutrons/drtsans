# Move part of the methods from api.py to avoid importing in loops
from mantid.simpleapi import mtd, logger, SaveAscii, RebinToWorkspace, SaveNexus  # noqa E402
# Import rolled up to complete a single top-level API
from drtsans import (apply_sensitivity_correction, solid_angle_correction)  # noqa E402
from drtsans import subtract_background  # noqa E402
from drtsans.transmission import apply_transmission_correction  # noqa E402
from drtsans.tof.eqsans.transmission import calculate_transmission  # noqa E402
from drtsans.thickness_normalization import normalize_by_thickness  # noqa E402
from drtsans.tof.eqsans.dark_current import subtract_dark_current  # noqa E402
from drtsans.mask_utils import apply_mask   # noqa E402
from drtsans.tof.eqsans.normalization import normalize_by_flux  # noqa E402
from drtsans.tof.eqsans.momentum_transfer import convert_to_q, split_by_frame  # noqa E402
from drtsans.dataobjects import IQmod, IQazimuthal
import os
import numpy as np
from typing import Tuple, List
from collections import namedtuple

# Binning parameters
BinningSetup = namedtuple('binning_setup', 'nxbins_main nybins_main n1dbins n1dbins_per_decade '
                                           'decade_on_center bin1d_type log_scale qmin, qmax, qxrange, qyrange')


# TODO - remove after 689 series is completely finished
def _convert_background_to_q(processed_background_ws, processed_sample_ws) -> Tuple[List[IQmod], List[IQazimuthal]]:
    """Re-process background to align binning with sample, convert to Q and split frames

    Parameters
    ----------
    processed_background_ws: mantid.dataobjects.EventWorkspace
        Processed background workspace
    processed_sample_ws: mantid.dataobjects.EventWorkspace
        Processed sample workspace

    Returns
    -------
    ~tuple
        list of I(Q1D), list of I(Q2D)

    """
    # Re-process background: workspace shall be rebinned to match the sample workspace
    bkgd_x0_vec = processed_background_ws.extractX()[0]
    sample_x0_vec = processed_sample_ws.extractX()[0]
    if len(bkgd_x0_vec) != len(sample_x0_vec):
        raise RuntimeError('Background workspace and sample workspace have different dimension on wavelength (X)')
    elif not np.allclose(sample_x0_vec, bkgd_x0_vec):
        # wavelength bins between sample and background are different
        # rebin background
        processed_background_ws = RebinToWorkspace(WorkspaceToRebin=processed_background_ws,
                                                   WorkspaceToMatch=processed_sample_ws,
                                                   OutputWorkspace=str(processed_background_ws))

    # convert to Q: Q1D and Q2D
    # No subpixel binning supported
    background_iq1d = convert_to_q(processed_background_ws, mode='scalar')
    background_iq2d = convert_to_q(processed_background_ws, mode='azimuthal')
    # split to frames
    background_iq1d_frames = split_by_frame(processed_background_ws, background_iq1d, verbose=True)
    background_iq2d_frames = split_by_frame(processed_background_ws, background_iq2d, verbose=True)

    return background_iq1d_frames, background_iq2d_frames


# # This is a composite method.  It can be used by
# # reduction_api.process_single_configuration_incoherence_correction()
# # without binning.
# def process_convert_q(raw_ws,
#                       transmission: Tuple[Any, float],
#                       theta_dependent_transmission,
#                       dark_current, flux, mask,
#                       solid_angle, sensitivity_workspace,
#                       sample_thickness: float,
#                       absolute_scale: float,
#                       output_suffix: str,
#                       delete_raw: bool) -> Tuple[List[IQmod], List[IQazimuthal], Any]:
#     """Process raw workspace and convert to Q and split into frames
#
#     Parameters
#     ----------
#     raw_ws:
#         raw event workspace and monitor workspace to process from
#     transmission: ~tuple
#         transmission workspace, transmission value
#     theta_dependent_transmission:
#         blabla
#     dark_current:
#         blabla
#     flux: ~tuple
#         flux method, flux run
#     mask: ~tuple
#         mask workspace, mask panel, mask BTP
#     solid_angle: bool
#         flag to do solid angle correction
#     sensitivity_workspace:
#         sensitivities workspace
#     sample_thickness: float
#         sample thickness in mm
#     absolute_scale: float
#         scale factor to intensities
#     output_suffix: float
#         suffix for output workspace
#     delete_raw: bool
#         flag to delete raw workspace
#
#     Returns
#     -------
#     ~tuple
#         list of IQmod, list of IQazimuthal, processed workspace
#
#     """
#     # Sanity check
#     assert raw_ws, 'Raw workspace cannot be None'
#
#     # Process raw workspace
#     output_workspace = str(raw_ws)
#     processed_ws = process_workspace_single_configuration(raw_ws, transmission, theta_dependent_transmission,
#                                                           dark_current, flux, mask,
#                                                           solid_angle, sensitivity_workspace,
#                                                           sample_thickness, absolute_scale,
#                                                           output_workspace, output_suffix)
#
#     # Optionally delete raw workspace
#     if delete_raw:
#         if isinstance(raw_ws, tuple):
#             raw_ws = raw_ws[0]
#         assert str(raw_ws) != str(processed_ws), 'Raw workspace and processed workspace have same name'
#         raw_ws.delete()
#
#     # No subpixel binning supported
#     # convert to Q: Q1D and Q2D
#     iq1d_main_in = convert_to_q(processed_ws, mode='scalar')
#     iq2d_main_in = convert_to_q(processed_ws, mode='azimuthal')
#     # split to frames
#     iq1d_main_in_fr = split_by_frame(processed_ws, iq1d_main_in, verbose=True)
#     iq2d_main_in_fr = split_by_frame(processed_ws, iq2d_main_in, verbose=True)
#
#     return iq1d_main_in_fr, iq2d_main_in_fr, processed_ws


# def process_workspace_single_configuration(ws_raw,
#                                            transmission,
#                                            theta_dependent_transmission,
#                                            dark_current,
#                                            flux,
#                                            mask,
#                                            solid_angle,
#                                            sensitivity_workspace,
#                                            thickness=1.,
#                                            absolute_scale=1.,
#                                            output_workspace=None,
#                                            output_suffix=''):
#     """This function provides quasi-full data processing for a single experimental configuration,
#     starting from workspaces (no data loading is happening inside this function)
#
#     This is a simplified version of eqsans.api.process_single_configuration().
#     The major difference is that
#     1. this method does not assume input workspace is a sample run
#     2. this method does not remove background
#     3. this method tends to use a more concise list of input parameters
#
#     Parameters
#     ----------
#     ws_raw: namedtuple
#         (~mantid.dataobjects.Workspace2D, ~mantid.dataobjects.Workspace2D)
#         raw data histogram workspace and monitor
#     sample_trans_ws:  ~mantid.dataobjects.Workspace2D
#         optional histogram workspace for sample transmission (already prepared)
#     sample_trans_value: float
#         optional value for sample transmission
#     bkg_ws_raw: ~mantid.dataobjects.Workspace2D
#         optional raw histogram workspace for background
#     bkg_trans_ws: ~mantid.dataobjects.Workspace2D
#         optional histogram workspace for background transmission
#     bkg_trans_value: float
#         optional value for background transmission
#     theta_deppendent_transmission: bool
#         flag to apply angle dependent transmission
#     dark_current: ~mantid.dataobjects.Workspace2D
#         dark current workspace
#     flux_method: str
#         normalization by time or monitor
#     mask_ws: ~mantid.dataobjects.Workspace2D
#         user defined mask
#     mask_panel: str
#         mask fron or back panel
#     mask_btp: dict
#         optional bank, tube, pixel to mask
#     solid_angle: bool
#         flag to apply solid angle
#     sensitivity_workspace: ~mantid.dataobjects.Workspace2D
#         workspace containing sensitivity
#     output_workspace: str
#         output workspace name
#     output_suffix:str
#         suffix for output workspace
#     thickness: float
#         sample thickness (cm)
#     absolute_scale_method: str
#         method to do absolute scaling (standard or direct_beam)
#     empty_beam_ws: ~mantid.dataobjects.Workspace2D
#         empty beam workspace for absolute scaling
#     beam_radius: float, None
#         beam radius for absolute scaling
#     absolute_scale: float
#         absolute scaling value for standard method
#     keep_processed_workspaces: bool
#         flag to keep the processed background workspace
#
#     Returns
#     -------
#     ~mantid.dataobjects.Workspace2D
#         Reference to the processed workspace
#     """
#     # Default output workspace name
#     if not output_workspace:
#         output_workspace = f'{str(ws_raw)}_single_config_{output_suffix}'
#
#     # Process input function parameters
#     flux_method, flux_value = flux
#     mask_ws, mask_panel, mask_btp = mask
#
#     # Prepare data workspace with dark current, flux, mask, solid angle and sensitivities
#     # create a common configuration for prepare data
#     prepare_data_conf = {'dark_current': dark_current,
#                          'flux_method': flux_method,
#                          'flux': flux_value,
#                          'mask_ws': mask_ws,
#                          'mask_panel': mask_panel,
#                          'mask_btp': mask_btp,
#                          'solid_angle': solid_angle,
#                          'sensitivity_workspace': sensitivity_workspace}
#     raw_ws = prepare_data_workspaces(ws_raw,
#                                      output_workspace=output_workspace,
#                                      **prepare_data_conf)
#
#     # Apply transmission to the sample
#     sample_trans_ws, sample_trans_value = transmission
#     print(f'tpe of transmission: {type(transmission)}')
#     if sample_trans_ws or sample_trans_value:
#         print(f'sample trans ws : {sample_trans_ws}\n\t\ttype = {type(sample_trans_ws)}')
#         print(f'sample trans val: {sample_trans_value}\n\t\ttype = {type(sample_trans_value)}')
#         if sample_trans_ws:
#             RebinToWorkspace(WorkspaceToRebin=sample_trans_ws,
#                              WorkspaceToMatch=raw_ws,
#                              OutputWorkspace=sample_trans_ws)
#         raw_ws = apply_transmission_correction(raw_ws,
#                                                trans_workspace=sample_trans_ws,
#                                                trans_value=sample_trans_value,
#                                                theta_dependent=theta_dependent_transmission,
#                                                output_workspace=output_workspace)
#
#     # finalize with absolute scale and thickness
#     # TODO FIXME @changwoo - normalization of sample thickness shall be applied to sample/background before or after
#     # inelastic/incoherence correction????
#     # Mathematically it does not mastter
#     raw_ws = normalize_by_thickness(raw_ws, thickness)
#     # absolute scale
#     raw_ws *= absolute_scale
#
#     return raw_ws


def prepare_data_workspaces(data: namedtuple,
                            dark_current=None,
                            flux_method=None,    # normalization (proton charge/time/monitor)
                            flux=None,           # additional file for normalization
                            mask_ws=None,        # apply a custom mask from workspace
                            mask_panel=None,     # mask back or front panel
                            mask_btp=None,       # mask bank/tube/pixel
                            solid_angle=True,
                            sensitivity_workspace=None,
                            output_workspace=None):

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

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = str(data.data)
        output_workspace = output_workspace.replace('_raw_histo', '') + '_processed_histo'

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


# NOTE: transformed from block of codes inside reduce_single_configuration
#       for calculating transmission
def process_transmission(transmission_ws, empty_trans_ws, transmission_radius, sensitivity_ws,
                         flux_method, flux, prefix,
                         type_name, output_dir, output_file_name):
    # sample transmission
    processed_transmission_dict = {}  # for output log
    raw_transmission_dict = {}  # for output log

    if transmission_ws.data is not None and empty_trans_ws is not None:
        # process transition workspace from raw
        processed_trans_ws_name = f'{prefix}_{type_name}_trans'  # type_name: sample/background
        processed_trans_ws = prepare_data_workspaces(transmission_ws,
                                                     flux_method=flux_method,
                                                     flux=flux,
                                                     solid_angle=False,
                                                     sensitivity_workspace=sensitivity_ws,
                                                     output_workspace=processed_trans_ws_name)
        # calculate transmission with fit function (default) Formula=a*x+b'
        calculated_trans_ws = calculate_transmission(processed_trans_ws, empty_trans_ws,
                                                     radius=transmission_radius, radius_unit="mm")
        print(f'{type_name} transmission =', calculated_trans_ws.extractY()[0, 0])

        # optionally save
        if output_dir:
            # save calculated transmission
            transmission_filename = os.path.join(output_dir, f'{output_file_name}_trans.txt')
            SaveAscii(calculated_trans_ws, Filename=transmission_filename)
            # Prepare result for drtsans.savereductionlog
            processed_transmission_dict['value'] = calculated_trans_ws.extractY()
            processed_transmission_dict['error'] = calculated_trans_ws.extractE()
            processed_transmission_dict['wavelengths'] = calculated_trans_ws.extractX()

            # Prepare result for drtsans.savereductionlog including raw sample transmission
            sample_trans_raw_ws = calculate_transmission(processed_trans_ws, empty_trans_ws,
                                                         radius=transmission_radius, radius_unit="mm",
                                                         fit_function='')

            raw_tr_fn = os.path.join(output_dir, f'{output_file_name}_raw_trans.txt')
            SaveAscii(sample_trans_raw_ws, Filename=raw_tr_fn)
            # Prepare result for drtsans.savereductionlog
            raw_transmission_dict['value'] = sample_trans_raw_ws.extractY()
            raw_transmission_dict['error'] = sample_trans_raw_ws.extractE()
            raw_transmission_dict['wavelengths'] = sample_trans_raw_ws.extractX()
    else:
        calculated_trans_ws = None

    return calculated_trans_ws, processed_transmission_dict, raw_transmission_dict
