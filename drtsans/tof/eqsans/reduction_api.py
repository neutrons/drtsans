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
from drtsans.tof.eqsans.correction_api import (CorrectionConfiguration, bin_i_of_q_per_wavelength,
                                               process_convert_q, do_inelastic_incoherence_correction_q1d)
import os
import numpy as np
from typing import Tuple, Any, List
from collections import namedtuple


# Binning parameters
binning_setup = namedtuple('binning_setup', 'nxbins_main nybins_main n1dbins n1dbins_per_decade '
                                            'decade_on_center bin1d_type log_scale qmin, qmax, qxrange, qyrange')


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
    print(f'[DEBUG WORKFLOW] Prove processed workspace: ')
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

    # Re-process background to align binning with sample, convert to Q and split frames
    processed_bkgd_ws = bkgd_raw_iq[2]
    background_iq1d_frames, background_iq2d_frames = _convert_background_to_q(processed_bkgd_ws,
                                                                              processed_sample_ws)

    # Subtract background from sample in workspace level for future output
    pure_sample_ws = subtract_background(processed_sample_ws, processed_bkgd_ws)

    # Process each frame separately!
    # output
    binned_iq1d_frames = list()
    binned_iq2d_frames = list()
    frame_q_range = list()

    replace_qmin = binning_params.qmin is None
    replace_qmax = binning_params.qmax is None
    n_wl_frames = len(raw_q1d_frame)

    print(f'[DEBUG INFO] Number of frames = {n_wl_frames}')

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

        # Bin sample data (without background), background and optionally elastic reference separately
        binned_sample_iq = bin_i_of_q_per_wavelength(raw_iq1d, raw_iq2d, binning_params)
        binned_bkgd_iq = bin_i_of_q_per_wavelength(background_iq1d_frames[frame], background_iq2d_frames[frame],
                                                   binning_params)

        if incoherence_correction_setup.elastic_reference_run:
            raw_ref_iq1d, raw_ref_iq2d = incoherence_correction_setup.elastic_reference_run.binned_ref_iq[frame]
            binned_elastic_ref_iq = bin_i_of_q_per_wavelength(raw_ref_iq1d, raw_ref_iq2d, binning_params)
        else:
            binned_elastic_ref_iq = None

        # step 4. process data with incoherent/inelastic correction
        if not incoherence_correction_setup.do_correction:
            # do not do any correction as a test for refactored reduction workflow, which shall
            # yield the same result if correction is skipped.
            binned_bkgd_iq1d, binned_bkgd_iq2d = binned_bkgd_iq
            binned_sample_iq1d, binned_sample_iq2d = binned_sample_iq
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

            print(f'[DEBUG 689] (Step 3 and 4) correct sample and background with {incoherence_correction_setup}')

            # correct 1D
            # correct I and dI of background accounting wavelength-dependent incoherent/inelastic scattering
            returned_list = do_inelastic_incoherence_correction_q1d([binned_sample_iq[0], binned_bkgd_iq[0]],
                                                                    incoherence_correction_setup)
            corrected_sample_1d, corrected_bkgd_1d = returned_list
            binned_sample_iq1d = corrected_sample_1d.iq1d
            binned_bkgd_iq1d = corrected_bkgd_1d.iq1d

            # Leave this to Jesse
            binned_sample_iq2d = binned_sample_iq[1]
            binned_bkgd_iq2d = binned_bkgd_iq[1]

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

        # Process 1D case
        # binned_sample_q1d = subtract_background(binned_sample_iq[0], binned_bkgd_iq[0])
        binned_sample_q1d = subtract_background(binned_sample_iq1d, binned_bkgd_iq1d)

        print(f'[DEBUG] sample q1D workspace: {binned_sample_q1d}')
        # TODO NOW - build IQmod from workspace
        vec_q = binned_sample_q1d.extractX().flatten()
        vec_i = binned_sample_q1d.extractY().flatten()
        vec_e = binned_sample_q1d.extractE().flatten()
        print(f'Q vec: size {vec_q.shape}')
        print(f'original IQ vec: size {binned_sample_iq[0].mod_q.shape}')
        # FIXME : q1d, wavelength order is different
        #         np.testing.assert_allclose(binned_sample_iq[0].mod_q, vec_q)
        print(f'NaN in I(Q): {len(np.where(np.isnan(vec_i))[0])}')
        print(f'NaN in I(Q): {len(np.where(np.isnan(binned_sample_iq[0].intensity))[0])}')

        # Separate wavelength, reconstruct to IQmod and remove nan and inf
        binned_sample_q1d = IQmod(vec_i, vec_e, vec_q, wavelength=binned_sample_iq[0].wavelength)
        binned_sample_q1d = binned_sample_q1d.be_finite()

        # Process 2D case
        binned_sample_q2d = subtract_background(binned_sample_iq2d, binned_bkgd_iq2d)
        # TODO FIXME - #743 binned_sample_q2d = binned_sample_q2d.be_finite()

        # append
        binned_iq1d_frames.append(binned_sample_q1d)
        binned_iq2d_frames.append(binned_sample_q2d)
    # END-FOR-FRAME

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


def process_elastic_reference_data(elastic_ref_setup, transmission_radius, sensitivity_ws, flux):
    """Process elastic reference run from raw workspaces to I of Q1D and Q2D split to frames

    Workflow
    1. Process ref sample transmission
    2. Process ref background transmission
    3. Process (single configuration) ref sample run
    4. Process (single configuration) ref background run
    5. Ref sample run convert to Q and split frame
    6. Ref background run convert to Q and split frame

    Requires:
      - dark current
      - mask
      - sensitivities
      - empty beam
      - empty beam radius

    Parameters
    ----------
    elastic_ref_setup: ElasticReferenceRunSetup
        elastic scattering correction reference setup

    Returns
    -------
    ElasticReferenceRunSetup


    """
    sample_transmission = elastic_ref_setup.ref_transmission_ws
    empty_trans_ws = elastic_ref_setup.empty_trans_ws

    sample_cal_trans_ws, _, _ = process_transmission(sample_transmission, empty_trans_ws,
                                                     transmission_radius, sensitivity_ws,
                                                     flux.method, flux.data, 'elastic_reference',
                                                     'elastic_reference', None, None)
    # process raw.
    assert sample_cal_trans_ws == 1, 'IMPLEMENTATION TO BE CONTINUED ... '

    # # process transmission if there is any
    # if sample_transmission.data is not None and empty_trans_ws is not None:
    #     sample_trans_ws_name = f'{prefix}_sample_trans'
    #     sample_trans_ws_processed = prepare_data_workspaces(loaded_ws.sample_transmission,
    #                                                         flux_method=flux_method,
    #                                                         flux=flux,
    #                                                         solid_angle=False,
    #                                                         sensitivity_workspace=loaded_ws.sensitivity,
    #                                                         output_workspace=sample_trans_ws_name)
    #     # calculate transmission with fit function (default) Formula=a*x+b'
    #     trans_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
    #                                       radius=transmission_radius, radius_unit="mm")

    return elastic_ref_setup
