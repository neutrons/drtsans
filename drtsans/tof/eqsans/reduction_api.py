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


def process_single_configuration_incoherence_correction(incoherence_correction_setup):
    """Process raw sample workspace with single configuration and inelastic/incoherence correction
    till binned I(Q, wavelength)

    Parameters
    ----------
    incoherence_correction_setup: CorrectionConfiguration
        Incoherence correction setup

    Returns
    -------

    """

    # 1. process single configuration of a sample run
    pass

    # 2. determine Qmin and Qmax sample run
    # min_q = max_q = None

    # 3. bin sample data (without background)
    # process data with incoherent/inelastic correction
    # sample_1d_fr, sample_2d_fr = process_bin_workspace(raw_sample_ws,
    #                                                    (sample_trans_ws, sample_trans_value),
    #                                                    theta_deppendent_transmission,
    #                                                    loaded_ws.dark_current,
    #                                                    (flux_method, flux),
    #                                                    (loaded_ws.mask, mask_panel, None),
    #                                                    solid_angle,
    #                                                    loaded_ws.sensitivity,
    #                                                    incoherence_correction_setup.sample_thickness,
    #                                                    absolute_scale,
    #                                                    binning_params)
    #
    # # FIXME - this if-else block is for debugging refined workflow.
    # if incoherence_correction_setup.debug_no_correction is False:
    #     # normalize
    #     if norm_dict:
    #         sample_1d_fr, sample_2d_fr = normalize_ws_with_elastic_scattering(sample_1d_fr,
    #                                                                           sample_2d_fr,
    #                                                                           norm_dict)
    #     # correct I and dI of background accounting wavelength-dependent incoherent/inelastic scattering
    #     r = correct_acc_incoherence_scattering(sample_1d_fr, sample_2d_fr, incoherence_correction_setup)
    #     iq1d_main_in_fr, iq2d_main_in_fr = r
    # else:
    #     # no correction for debugging purpose
    #     # FIXME - remove this after everything passes!
    #     iq1d_main_in_fr = sample_1d_fr
    #     iq2d_main_in_fr = sample_2d_fr
    #
    # # subtract with background
    # print(f'Binning: {binning_params}')
    # print(f'Number of frames: {len(iq1d_main_in_fr)}')
    # for i_f in range(len(iq1d_main_in_fr)):
    #     print('1D')
    #     print(f'[NOW-CORRECTION] 1D: sample     range {iq1d_main_in_fr[i_f].mod_q[0]}, '
    #           f'{iq1d_main_in_fr[i_f].mod_q[-1]}')
    #     print(f'[NOW-CORRECTION] 1D: background range {bkgd_iq1d[i_f].mod_q[0]}, '
    #           f'{bkgd_iq1d[i_f].mod_q[-1]}')
    #     iq1d_main_in_fr[i_f] = subtract_background(iq1d_main_in_fr[i_f], bkgd_iq1d[i_f])
    #
    #     print('2D')
    #     print(f'[NOW-CORRECTION] 2D: range {iq2d_main_in_fr.qx[0, 0]}, '
    #           f'{iq2d_main_in_fr.qx[0, nxbins_main - 1]}')
    #     iq2d_main_in_fr[i_f] = subtract_background(iq2d_main_in_fr[i_f], bkgd_iq2d[i_f])
    #     iq2d_main_in_fr[i] = iq2d_main_in_fr[i] - bkgd_iq2d[i]

    # 4. bin background data
    pass

    # 5. normalize by elastic scattering reference data
    if incoherence_correction_setup.do_correction and incoherence_correction_setup.elastic_reference_run:
        # normalize
        # calculate the normalization factors
        # # calculate normalization factor
        # norm_dict = calculate_elastic_scattering_factor(loaded_ws.elastic_ref_ws,
        #                                                 loaded_ws.elastic_ref_trans_ws,
        #                                                 elastic_ref_setup.ref_trans_value,
        #                                                 elastic_ref_setup.ref_sample_thickness,
        #                                                 None)

        # normalize background and sample
        # bkgd_iq1d, bkgd_iq2d = normalize_ws_with_elastic_scattering(bkgd_iq1d, bkgd_iq2d, norm_dict)
        # sample_iq1d, sample_iq2d = normalize_ws_with_elastic_scattering(sample_iq1d, sample_iq2d, norm_dict)
        print(f'ASAP')

    # 6. correct to reduce incoherence and inelastic scattering
    if incoherence_correction_setup.do_correction:
        # bkgd_iq1d, bkgd_iq2d = correct_acc_incoherence_scattering(bkgd_iq1d, bkgd_iq2d,
        #                                                           incoherence_correction_setup)
        # sample_iq1d, sample_iq2d = correct_acc_incoherence_scattering(bkgd_iq1d, bkgd_iq2d,
        #                                                           incoherence_correction_setup)
        pass

    # 7. Remove background
    # subtract_bkgd()

    # 8. Return Q1D and Q2D

    return None, None


# TODO - compare this method with process_raw_workspace()
def process_background_workspace():
    """Process raw background workspace including
    (1) regular process data workspace
    (2) apply transmission correction
    (3) bin to Q
    (4) split to frame

    Product:
    1. returned Q frames
    2. processed background workspace

    Returns
    -------

    """
    # # Determine workspace name
    # bkgd_ws_name = output_suffix + '_background'
    # bkgd_ws = prepare_data_workspaces(bkg_ws_raw,
    #                                   output_workspace=bkgd_ws_name,
    #                                   **prepare_data_conf)
    # # apply transmission to bkgd
    # if bkg_trans_ws or bkg_trans_value:
    #     # make binning consistent
    #     if bkg_trans_ws:
    #         RebinToWorkspace(WorkspaceToRebin=bkg_trans_ws,
    #                          WorkspaceToMatch=bkgd_ws,
    #                          OutputWorkspace=bkg_trans_ws)
    #     # apply transmission correction to output workspace
    #     bkgd_ws = apply_transmission_correction(bkgd_ws,
    #                                             trans_workspace=bkg_trans_ws,
    #                                             trans_value=bkg_trans_value,
    #                                             theta_dependent=theta_dependent_transmission,
    #                                             output_workspace=bkgd_ws_name)
    # # normalize by thickness
    # bkgd_ws = normalize_by_thickness(bkgd_ws, thickness)
    #
    # # scale up
    # bkgd_ws *= absolute_scale
    #
    # # convert to Q
    # iq1d_main_in = convert_to_q(bkgd_ws, mode='scalar')
    # iq2d_main_in = convert_to_q(bkgd_ws, mode='azimuthal')
    #
    # # split to frames
    # iq1d_main_in_fr = split_by_frame(bkgd_ws, iq1d_main_in)
    # iq2d_main_in_fr = split_by_frame(bkgd_ws, iq2d_main_in)

    iq1d_main_in_fr = iq2d_main_in_fr = None

    return iq1d_main_in_fr, iq2d_main_in_fr


def process_raw_workspace(ws_raw,
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
    if not output_workspace:
        output_workspace = output_suffix + '_sample'

    # expand
    flux_method, flux_value = flux
    mask_ws, mask_panel, mask_btp = mask

    # create a common configuration for prepare data
    prepare_data_conf = {'dark_current': dark_current,
                         'flux_method': flux_method,
                         'flux': flux_value,
                         'mask_ws': mask_ws,
                         'mask_panel': mask_panel,
                         'mask_btp': mask_btp,
                         'solid_angle': solid_angle,
                         'sensitivity_workspace': sensitivity_workspace}

    # process sample
    raw_ws = prepare_data_workspaces(ws_raw,
                                     output_workspace=output_workspace,
                                     **prepare_data_conf)
    # apply transmission to the sample
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
    raw_ws = normalize_by_thickness(raw_ws, thickness)
    # absolute scale
    raw_ws *= absolute_scale

    return raw_ws


# TODO/FIXME - this is an exact copy of api.prepare_data_workspaces
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
