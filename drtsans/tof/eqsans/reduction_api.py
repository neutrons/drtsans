# Move part of the methods from api.py to avoid importing in loops
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
from drtsans.dataobjects import IQmod, IQazimuthal
from drtsans.tof.eqsans.correction_api import bin_i_of_q
import os
from drtsans.tof.eqsans.correction_api import process_convert_q
import numpy as np
from typing import Tuple, Any, List
from collections import namedtuple


# Binning parameters
binning_setup = namedtuple('binning_setup', 'nxbins_main nybins_main n1dbins n1dbins_per_decade decade_on_center bin1d_type log_scale qmin, qmax, qxrange, qyrange')


def process_single_configuration_incoherence_correction(sample_ws, sample_transmission,
                                                        theta_dependent_transmission,
                                                        dark_current,
                                                        flux_setup,
                                                        mask_setup,
                                                        solid_angle,
                                                        sensitivities_ws,
                                                        absolute_scale,
                                                        sample_thickness,
                                                        bkgd_raw_iq: Tuple[List[IQmod], List[IQazimuthal], Any],
                                                        incoherence_correction_setup,
                                                        binning_params: binning_setup
                                                        ) -> Tuple[List[Any], List[Any], Any,
                                                                   List[binning_setup]]:
    """Process raw sample workspace with single configuration and inelastic/incoherence correction
    till binned I(Q, wavelength)

    Parameters
    ----------
    bkgd_raw_iq: ~tuple
        List of I(Q), List of I(Qx, Qy) of various frames, Processed background workspace
    incoherence_correction_setup: CorrectionConfiguration
        Incoherence correction setup
    binning_params: binning_setup
        binning parameters including range, step and type

    Returns
    -------
    tuple
        list (binned Q1D), list (binned Q2D), processed workspace, list of tuples as range of Q, Qx, Qy

    """
    # Check input parameters
    assert isinstance(incoherence_correction_setup, CorrectionConfiguration)

    # Process single configuration of a sample run (raw)
    sample_raw_iq = process_convert_q(sample_ws,
                                      sample_transmission,  # (sample_trans_ws, sample_trans_value),
                                      theta_dependent_transmission,
                                      dark_current,
                                      flux_setup,  # (flux_method, flux),
                                      mask_setup,  # (loaded_ws.mask, mask_panel, None),
                                      solid_angle,
                                      sensitivities_ws,
                                      sample_thickness,
                                      absolute_scale,
                                      'sample',
                                      delete_raw=True)
    raw_q1d_frame, raw_q2d_frame, processed_sample_ws = sample_raw_iq

    # Re-process background: workspace shall be rebinned to match the sample workspace
    print(f'[DEBUG WORKFLOW] Prove processed workspace: ')
    processed_bkgd_ws = bkgd_raw_iq[2]
    bkgd_x0_vec = processed_bkgd_ws.extractX()[0]
    sample_x0_vec = processed_sample_ws.extractX()[0]
    if len(bkgd_x0_vec) != len(sample_x0_vec):
        raise RuntimeError('Background workspace and sample workspace have different dimension on wavelength (X)')
    elif not np.allclose(sample_x0_vec, bkgd_x0_vec):
        # wavelength bins between sample and background are different
        # rebin background
        processed_bkgd_ws = RebinToWorkspace(WorkspaceToRebin=processed_bkgd_ws,
                                             WorkspaceToMatch=processed_sample_ws,
                                             OutputWorkspace=str(processed_bkgd_ws))

    # convert to Q: Q1D and Q2D
    # No subpixel binning supported
    background_iq1d = convert_to_q(processed_bkgd_ws, mode='scalar')
    background_iq2d = convert_to_q(processed_bkgd_ws, mode='azimuthal')
    # split to frames
    background_iq1d_frames = split_by_frame(processed_bkgd_ws, background_iq1d, verbose=True)
    background_iq2d_frames = split_by_frame(processed_bkgd_ws, background_iq2d, verbose=True)

    # FIXME - DEBUG OUTPUT
    print(f'sample     workspace unit X: {processed_sample_ws.getAxis(0).getUnit().unitID()}')
    print(f'background workspace unit X: {processed_bkgd_ws.getAxis(0).getUnit().unitID()}')
    wl_array = processed_sample_ws.extractX()
    print(f'sample     wavelength size: {wl_array.shape}')
    wl_array = processed_bkgd_ws.extractX()
    print(f'background wavelength size: {wl_array.shape}')
    wl_array = processed_sample_ws.extractX()
    print(f'sample     wavelength range: '
          f'{wl_array[0][0]}, {wl_array[0][-1]} ... {wl_array[-1][0]}, {wl_array[-1][-1]}')
    wl_array = processed_bkgd_ws.extractX()
    print(f'background wavelength range: '
          f'{wl_array[0][0]}, {wl_array[0][-1]} ... {wl_array[-1][0]}, {wl_array[-1][-1]}')
    # END ====================================================================================

    # Subtract background from sample in workspace level for future output
    pure_sample_ws = subtract_background(processed_sample_ws, processed_bkgd_ws)
    debug_iq1d = convert_to_q(pure_sample_ws, mode='scalar')
    print(f'[DEBUG] Processed sample - background.  Q range = {debug_iq1d.mod_q.min()}, {debug_iq1d.mod_q.max()}')

    # # convert to Q
    # iq1d_main_in = convert_to_q(pure_sample_ws, mode='scalar')
    # iq2d_main_in = convert_to_q(pure_sample_ws, mode='azimuthal')
    # # split to frames
    # iq1d_main_in_fr = split_by_frame(pure_sample_ws, iq1d_main_in)
    # iq2d_main_in_fr = split_by_frame(pure_sample_ws, iq2d_main_in)
    # get frames
    # n_wl_frames = len(iq2d_main_in_fr)
    # # Process each frame separately

    # # Step 2 to 4 are up to each frame
    # # Determine Q range for each frame from PURE sample I(Q1D) and I(Q2D)
    # frame_q_range = list()
    # for wl_frame in range(n_wl_frames):
    #     iq2d = iq2d_main_in_fr[wl_frame]
    #     iq1d = iq1d_main_in_fr[wl_frame]
    #     mod_q_range = iq1d.mod_q.min(), iq1d.mod_q.max()
    #     qx_range = iq2d.qx.min(), iq2d.qx.max()
    #     qy_range = iq2d.qy.min(), iq2d.qy.max()
    #     frame_q_range.append((mod_q_range, qx_range, qy_range))
    #     print(f'[DEBUG PROCESSED FRAME] Frame = {wl_frame}: Q1D range = {iq1d.mod_q.min()},'
    #           f'{iq1d.mod_q.max()}; Q2D X range = {iq2d.qx.min(), iq2d.qx.max()}, Y range = '
    #           f'{iq2d.qy.min()}, {iq2d.qy.max()},  Number data points (1D) = {len(iq1d.mod_q)}')
    # #

    # Process each frame separately!
    # output
    binned_iq1d_frames = list()
    binned_iq2d_frames = list()
    frame_q_range = list()

    replace_qmin = binning_params.qmin is None
    replace_qmax = binning_params.qmax is None
    n_wl_frames = len(raw_q1d_frame)
    for frame in range(n_wl_frames):
        # Bin raw I(Q), I(Qx, Qy), elastic normalization and inelastic/incoherent correction
        # step 2. determine binning parameters:  Qmin and Qmax from sample run
        # 1D
        raw_iq1d = raw_q1d_frame[frame]
        assert isinstance(raw_iq1d, IQmod), f'Raw I(Q1D) is of type {type(raw_iq1d)}'

        # Set Qmin and Qmax if they are not implemented
        if replace_qmin:
            binning_params = binning_params._replace(qmin=raw_iq1d.mod_q.min())
        if replace_qmax:
            binning_params = binning_params._replace(qmax=raw_iq1d.mod_q.max())

        # 2D
        raw_iq2d = raw_q2d_frame[frame]
        assert isinstance(raw_iq2d, IQazimuthal), f'Raw I(Q2D) is of type {type(raw_iq2d)}'

        if binning_params.qxrange is None:
            # default: data's qx range
            qx_min = np.min(raw_iq2d.qx)
            qx_max = np.max(raw_iq2d.qx)
            binning_params = binning_params._replace(qxrange=(qx_min, qx_max))

        if binning_params.qyrange is None:
            # default: data's qy range
            qy_min = np.min(raw_iq2d.qy)
            qy_max = np.max(raw_iq2d.qy)
            binning_params = binning_params._replace(qyrange=(qy_min, qy_max))

        # Set to future record
        frame_q_range.append(binning_params)

        print(f'[SAMPLE FRAME] {frame} Q range: {raw_iq1d.mod_q.min()}, {raw_iq1d.mod_q.max()}.  Number data points = {len(raw_iq1d.mod_q)}')
        print(f'[DEBUG binning] {binning_params.qxrange}, {binning_params.qyrange}')

        # Bin sample data (without background), background and optionally elastic reference separately
        binned_sample_iq = bin_i_of_q(raw_iq1d, raw_iq2d, binning_params)
        binned_bkgd_iq = bin_i_of_q(background_iq1d_frames[frame], background_iq2d_frames[frame], binning_params)

        if incoherence_correction_setup.elastic_reference_run:
            raw_ref_iq1d, raw_ref_iq2d = incoherence_correction_setup.elastic_reference_run.binned_ref_iq[frame]
            binned_elastic_ref_iq = bin_i_of_q(raw_ref_iq1d, raw_ref_iq2d, binning_params)
        else:
            binned_elastic_ref_iq = None

        # step 4. process data with incoherent/inelastic correction
        # FIXME - this if-else block is for debugging refined workflow.
        if incoherence_correction_setup.debug_no_correction:
            # do not do any correction as a test for refactored reduction workflow, which shall
            # yield the same result if correction is skipped.
            pass
        else:
            # normalize sample run and background run
            if binned_elastic_ref_iq:
                # calculate normalization according to new binning
                raise NotImplementedError('ASAP')
                # norm_dict = calculate_elastic_scattering_factor(loaded_ws.elastic_ref_ws,
                #                                                 loaded_ws.elastic_ref_trans_ws,
                #                                                 elastic_ref_setup.ref_trans_value,
                #                                                 elastic_ref_setup.ref_sample_thickness,
                #                                                 None)

                # normalize sample data
                # sample_iq1d, sample_iq2d = normalize_ws_with_elastic_scattering(sample_iq1d, sample_iq2d, norm_dict)

                # normalize background data
                # bkgd_iq1d, bkgd_iq2d = normalize_ws_with_elastic_scattering(bkgd_iq1d, bkgd_iq2d, norm_dict)

                # sample_1d_fr, sample_2d_fr = normalize_ws_with_elastic_scattering(sample_1d_fr,
                #                                                               sample_2d_fr,
                #                                                               norm_dict)
            # END-IF (correction step 2)

            # correct 1D
            # correct I and dI of background accounting wavelength-dependent incoherent/inelastic scattering
            # bkgd_iq1d, bkgd_iq2d = correct_acc_incoherence_scattering(bkgd_iq1d, bkgd_iq2d,
            #                                                           incoherence_correction_setup)
            # sample_iq1d, sample_iq2d = correct_acc_incoherence_scattering(bkgd_iq1d, bkgd_iq2d,
            #                                                           incoherence_correction_setup)
        # END-IF-step 4

        # Merge sample and background: subtract sample with background
        print(f'Binning: {binning_params}')
        print('1D {}, {}'.format(type(binned_sample_iq), type(binned_bkgd_iq)))
        print(f'[NOW-CORRECTION] 1D: sample     range {binned_sample_iq[0].mod_q[0]}, '
              f'{binned_sample_iq[0].mod_q[-1]}')
        print(f'[NOW-CORRECTION] 1D: background range {binned_bkgd_iq[0].mod_q[0]}, '
              f'{binned_bkgd_iq[0].mod_q[-1]}')
        print('2D')
        print(f'[NOW-CORRECTION] 2D: range {binned_bkgd_iq[1].qx[0, 0]}')

        binned_sample_q1d = subtract_background(binned_sample_iq[0], binned_bkgd_iq[0])
        print(f'[DEBUG] sample q1D workspace: {binned_sample_q1d}')
        # TODO NOW - build IQmod from workspace
        vec_q = binned_sample_q1d.extractX().flatten()
        vec_i = binned_sample_q1d.extractY().flatten()
        vec_e = binned_sample_q1d.extractE().flatten()
        print(f'Q vec: size {vec_q.shape}')
        print(f'original IQ vec: size {binned_sample_iq[0].mod_q.shape}')
        np.testing.assert_allclose(binned_sample_iq[0].mod_q, vec_q)
        print(f'NaN in I(Q): {len(np.where(np.isnan(vec_i))[0])}')
        print(f'NaN in I(Q): {len(np.where(np.isnan(binned_sample_iq[0].intensity))[0])}')

        # Separate wavelength, reconstruct to IQmod and remove nan and inf
        binned_sample_q1d = IQmod(vec_i, vec_e, vec_q, wavelength=binned_sample_iq[0].wavelength)
        binned_sample_q1d = binned_sample_q1d.be_finite()

        # 2D case
        binned_sample_q2d = subtract_background(binned_sample_iq[1], binned_bkgd_iq[1])
        try:
            binned_sample_q2d = binned_sample_q2d.be_finite()
        except AttributeError:
            # TODO FIXME - implement be_finite() for 2D
            pass

        # append
        binned_iq1d_frames.append(binned_sample_q1d)
        binned_iq2d_frames.append(binned_sample_q2d)

    return binned_iq1d_frames, binned_iq2d_frames, pure_sample_ws, frame_q_range


def process_workspace_single_configuration(ws_raw,
                                           transmission,
                                           theta_dependent_transmission,
                                           dark_current,
                                           flux,
                                           mask,
                                           solid_angle,
                                           sensitivity_workspace,
                                           thickness=1.,
                                           absolute_scale=1.,
                                           output_workspace=None,
                                           output_suffix=''):
    """This function provides quasi-full data processing for a single experimental configuration,
    starting from workspaces (no data loading is happening inside this function)

    This is a simplified version of eqsans.api.process_single_configuration().
    The major difference is that
    1. this method does not assume input workspace is a sample run
    2. this method does not remove background
    3. this method tends to use a more concise list of input parameters

    Parameters
    ----------
    ws_raw: namedtuple
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
    beam_radius: float, None
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
    # Default output workspace name
    if not output_workspace:
        output_workspace = f'{str(ws_raw)}_single_config_{output_suffix}'

    # Process input function parameters
    flux_method, flux_value = flux
    mask_ws, mask_panel, mask_btp = mask

    # Prepare data workspace with dark current, flux, mask, solid angle and sensitivities
    # create a common configuration for prepare data
    prepare_data_conf = {'dark_current': dark_current,
                         'flux_method': flux_method,
                         'flux': flux_value,
                         'mask_ws': mask_ws,
                         'mask_panel': mask_panel,
                         'mask_btp': mask_btp,
                         'solid_angle': solid_angle,
                         'sensitivity_workspace': sensitivity_workspace}
    raw_ws = prepare_data_workspaces(ws_raw,
                                     output_workspace=output_workspace,
                                     **prepare_data_conf)

    # Apply transmission to the sample
    sample_trans_ws, sample_trans_value = transmission
    print(f'tpe of transmission: {type(transmission)}')
    if sample_trans_ws or sample_trans_value:
        print(f'sample trans ws : {sample_trans_ws}\n\t\ttype = {type(sample_trans_ws)}')
        print(f'sample trans val: {sample_trans_value}\n\t\ttype = {type(sample_trans_value)}')
        if sample_trans_ws:
            RebinToWorkspace(WorkspaceToRebin=sample_trans_ws,
                             WorkspaceToMatch=raw_ws,
                             OutputWorkspace=sample_trans_ws)
        raw_ws = apply_transmission_correction(raw_ws,
                                               trans_workspace=sample_trans_ws,
                                               trans_value=sample_trans_value,
                                               theta_dependent=theta_dependent_transmission,
                                               output_workspace=output_workspace)

    # finalize with absolute scale and thickness
    # TODO FIXME @changwoo - normalization of sample thickness shall be applied to sample/background before or after
    # inelastic/incoherence correction????
    # Mathematically it does not mastter
    raw_ws = normalize_by_thickness(raw_ws, thickness)
    # absolute scale
    raw_ws *= absolute_scale

    return raw_ws


def prepare_data_workspaces(data,
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
