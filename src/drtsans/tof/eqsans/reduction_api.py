# Move part of the methods from api.py to avoid importing in loops
import os
from collections import namedtuple
from typing import Dict, List, Tuple

import numpy as np
from mantid.simpleapi import (
    SaveAscii,
    mtd,
)  # noqa E402

# Import rolled up to complete a single top-level API
from drtsans import (  # noqa E402
    apply_sensitivity_correction,
    solid_angle_correction,
    subtract_background,  # noqa E402
)
from drtsans.dataobjects import IQmod, IQazimuthal  # noqa E402
from drtsans.iq import bin_all  # noqa E402
from drtsans.mask_utils import apply_mask  # noqa E402
from drtsans.thickness_normalization import normalize_by_thickness  # noqa E402
from drtsans.tof.eqsans.correction_api import (
    CorrectionConfiguration,
    do_inelastic_incoherence_correction,
    save_k_vector,
)
from drtsans.tof.eqsans.dark_current import subtract_dark_current  # noqa E402
from drtsans.tof.eqsans.elastic_reference_normalization import (
    normalize_by_elastic_reference_all,
)
from drtsans.tof.eqsans.normalization import normalize_by_flux  # noqa E402
from drtsans.tof.eqsans.transmission import calculate_transmission  # noqa E402
from drtsans.transmission import apply_transmission_correction  # noqa E402

# Binning parameters
BinningSetup = namedtuple(
    "binning_setup",
    "nxbins_main nybins_main n1dbins n1dbins_per_decade "
    "decade_on_center bin1d_type log_scale qmin, qmax, qxrange, qyrange",
)


def prepare_data_workspaces(
    data: namedtuple,
    dark_current=None,
    flux_method=None,  # normalization (proton charge/time/monitor)
    flux=None,  # additional file for normalization
    mask_ws=None,  # apply a custom mask from workspace
    mask_panel=None,  # mask back or front panel
    mask_btp=None,  # mask bank/tube/pixel
    solid_angle=True,
    sensitivity_workspace=None,
    output_workspace=None,
    has_blocked_beam=False,
):
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
        output_workspace = output_workspace.replace("_raw_histo", "") + "_processed_histo"

    mtd[str(data.data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    # Dark current
    if dark_current is not None and dark_current.data is not None:
        subtract_dark_current(output_workspace, dark_current.data)

    # Normalization
    if flux_method is not None:
        kw = dict(method=flux_method, output_workspace=output_workspace, has_blocked_beam=has_blocked_beam)
        if flux_method == "monitor":
            kw["monitor_workspace"] = data.monitor
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
        apply_sensitivity_correction(output_workspace, sensitivity_workspace=sensitivity_workspace)

    return mtd[output_workspace]


# NOTE: transformed from block of codes inside reduce_single_configuration
#       for calculating transmission
def process_transmission(
    transmission_ws,
    empty_trans_ws,
    transmission_radius,
    sensitivity_ws,
    flux_method,
    flux,
    prefix,
    type_name,
    output_dir,
    output_file_name,
):
    # sample transmission
    processed_transmission_dict = {}  # for output log
    raw_transmission_dict = {}  # for output log

    if transmission_ws.data is not None and empty_trans_ws is not None:
        # process transition workspace from raw
        processed_trans_ws_name = f"{prefix}_{type_name}_trans"  # type_name: sample/background
        processed_trans_ws = prepare_data_workspaces(
            transmission_ws,
            flux_method=flux_method,
            flux=flux,
            solid_angle=False,
            sensitivity_workspace=sensitivity_ws,
            output_workspace=processed_trans_ws_name,
        )
        # calculate transmission with fit function (default) Formula=a*x+b'
        calculated_trans_ws = calculate_transmission(
            processed_trans_ws,
            empty_trans_ws,
            radius=transmission_radius,
            radius_unit="mm",
        )
        print(f"{type_name} transmission =", calculated_trans_ws.extractY()[0, 0])

        # optionally save
        if output_dir:
            # save calculated transmission
            transmission_filename = os.path.join(output_dir, f"{output_file_name}_trans.txt")
            SaveAscii(calculated_trans_ws, Filename=transmission_filename)
            # Prepare result for drtsans.savereductionlog
            processed_transmission_dict["value"] = calculated_trans_ws.extractY()
            processed_transmission_dict["error"] = calculated_trans_ws.extractE()
            processed_transmission_dict["wavelengths"] = calculated_trans_ws.extractX()

            # Prepare result for drtsans.savereductionlog including raw sample transmission
            sample_trans_raw_ws = calculate_transmission(
                processed_trans_ws,
                empty_trans_ws,
                radius=transmission_radius,
                radius_unit="mm",
                fit_function="",
            )

            raw_tr_fn = os.path.join(output_dir, f"{output_file_name}_raw_trans.txt")
            SaveAscii(sample_trans_raw_ws, Filename=raw_tr_fn)
            # Prepare result for drtsans.savereductionlog
            raw_transmission_dict["value"] = sample_trans_raw_ws.extractY()
            raw_transmission_dict["error"] = sample_trans_raw_ws.extractE()
            raw_transmission_dict["wavelengths"] = sample_trans_raw_ws.extractX()
    else:
        calculated_trans_ws = None

    return calculated_trans_ws, processed_transmission_dict, raw_transmission_dict


def bin_i_with_correction(
    iq1d_in_frames: List[IQmod],
    iq2d_in_frames: List[IQazimuthal],
    frameskip_frame: int,
    slice_name: str,
    weighted_errors: bool,
    user_qmin: float,
    user_qmax: float,
    num_x_bins: int,
    num_y_bins: int,
    num_q1d_bins: int,
    num_q1d_bins_per_decade: int,
    decade_on_center: bool,
    bin1d_type: str,
    log_binning: bool,
    annular_bin: float,
    wedges: List[Tuple[int, int]],
    symmetric_wedges: bool,
    correction_setup: CorrectionConfiguration,
    iq1d_elastic_ref_fr: List[IQmod],
    iq2d_elastic_ref_fr: List[IQazimuthal],
    raw_name: str,
    output_dir: str,
    output_filename: str = "",
) -> Tuple[IQazimuthal, List[IQmod]]:
    """Bin I(Q) in 1D and 2D with the option to do elastic and/or inelastic incoherent correction

    Parameters
    ----------
    iq1d_in_frames: list[~drtsans.dataobjects.IQmod]
        Objects containing 1D unbinned data I(\|Q\|). It will be used for scalar binned data
    iq2d_in_frames: list[~drtsans.dataobjects.IQazimuthal]
        Objects containing 2D unbinned data I(Qx, Qy). It will be used for 2D binned data,
        and 1D wedge or annular binned data
    frameskip_frame: int
        Index of the frame in ``iq1d_in_frames`` and ``iq2d_in_frames`` to bin and correct
    weighted_errors: bool
        If True, the binning is done using the Weighted method
    user_qmin: float
        Minimum value of the momentum transfer modulus Q
    user_qmax: float
        Maximum value of the momentum transfer modulus Q
    num_x_bins: int
        Number of bins in the x direction for 2D binning
    num_y_bins: int
        Number of bins in the y direction for 2D binning
    num_q1d_bins: int
        Number of bins for the 1d binning.
    num_q1d_bins_per_decade: int
        Total number of bins will be this value multiplied by
        number of decades from X min to X max
    decade_on_center: bool
        Flag to have the min X and max X on bin center; Otherwise, they will be on bin boundary
    bin1d_type: str
        Type of binning for 1D data. Possible choices are 'scalar', 'annular', or 'wedge'
    log_binning: bool
        If True, 1D scalar or wedge binning will be logarithmic. Ignored for anything else
    annular_bin: float
        Width of annular bin in degrees. Annular binning is linear
    wedges: list
        List of tuples (angle_min, angle_max) for the wedges. Both numbers have to be in
        the [-90,270) range. It will add the wedge offset by 180 degrees dependent
        on ``symmetric_wedges``
    symmetric_wedges: bool
        It will add the wedge offset by 180 degrees if True
    correction_setup: ~drtsans.tof.eqsans.correction_api.CorrectionConfiguration
        Parameters for elastic and inelastic/incoherence scattering correction
    iq1d_elastic_ref_fr: list[~drtsans.dataobjects.IQmod]
        Objects containing 1D unbinned data I(\|Q\|) for the elastic reference run
    iq2d_elastic_ref_fr: list[~drtsans.dataobjects.IQazimuthal]
        Objects containing 2D unbinned data I(Qx, Qy) for the elastic reference run
    raw_name: str
        Prefix for file to save inelastic incoherent correction (B) data
    output_dir: str
        Output directory for I(Q) profiles
    output_filename: str
        Output filename parsed from input configuration file (JSON)

    Returns
    -------
    (~drtsans.dataobjects.IQazimuthal, list[~drtsans.dataobjects.IQmod])
        - Binned and optionally corrected IQazimuthal
        - List of binned and optionally corrected IQmod objects. The list has length
          1, unless the 'wedge' mode is selected, when the length is the number of
          original wedges
    """

    # Setup for corrections
    if correction_setup.do_elastic_correction or any(correction_setup.do_inelastic_correction):
        # If any correction is turned on, then weighted_errors is always true
        weighted_errors = True

        # Define qmin and qmax for this frame
        if user_qmin is None:
            qmin = iq1d_in_frames[frameskip_frame].mod_q.min()
        else:
            qmin = user_qmin
        if user_qmax is None:
            qmax = iq1d_in_frames[frameskip_frame].mod_q.max()
        else:
            qmax = user_qmax

        # Set qxrange and qyrange for this frame
        qxrange = np.min(iq2d_in_frames[frameskip_frame].qx), np.max(iq2d_in_frames[frameskip_frame].qx)
        qyrange = np.min(iq2d_in_frames[frameskip_frame].qy), np.max(iq2d_in_frames[frameskip_frame].qy)

        # Bin I(Q1D, wl) and I(Q2D, wl) in Q and (Qx, Qy) space respectively but not wavelength
        iq2d_main_wl, iq1d_main_wl = bin_all(
            iq2d_in_frames[frameskip_frame],
            iq1d_in_frames[frameskip_frame],
            num_x_bins,
            num_y_bins,
            n1dbins=num_q1d_bins,
            n1dbins_per_decade=num_q1d_bins_per_decade,
            decade_on_center=decade_on_center,
            # corrections should use all the detector-panel's area, not just one wedge
            bin1d_type="scalar" if bin1d_type == "wedge" else bin1d_type,
            log_scale=log_binning,
            qmin=qmin,
            qmax=qmax,
            qxrange=qxrange,
            qyrange=qyrange,
            annular_angle_bin=annular_bin,
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            error_weighted=weighted_errors,
            n_wavelength_bin=None,
        )
        # Check due to functional limitation
        assert isinstance(iq1d_main_wl, list), f"Output I(Q) must be a list but not a {type(iq1d_main_wl)}"
        if len(iq1d_main_wl) != 1:
            raise NotImplementedError(f"Not expected that there are more than 1 IQmod main but {len(iq1d_main_wl)}")

    # Elastic correction
    if correction_setup.do_elastic_correction:
        elastic_output_dir = os.path.join(
            output_dir, "info", "elastic_norm", output_filename, slice_name, f"frame_{frameskip_frame}"
        )
        os.makedirs(elastic_output_dir, exist_ok=True)

        k_file_prefix = f"{raw_name}"

        # Bin elastic reference run
        if iq1d_elastic_ref_fr:
            # bin the reference elastic runs of the current frame
            iq2d_elastic_wl, iq1d_elastic_wl = bin_all(
                iq2d_elastic_ref_fr[frameskip_frame],
                iq1d_elastic_ref_fr[frameskip_frame],
                num_x_bins,
                num_y_bins,
                n1dbins=num_q1d_bins,
                n1dbins_per_decade=num_q1d_bins_per_decade,
                decade_on_center=decade_on_center,
                # corrections should use all the detector-panel's area, not just one wedge
                bin1d_type="scalar" if bin1d_type == "wedge" else bin1d_type,
                log_scale=log_binning,
                qmin=qmin,
                qmax=qmax,
                qxrange=qxrange,
                qyrange=qyrange,
                annular_angle_bin=annular_bin,
                wedges=wedges,
                symmetric_wedges=symmetric_wedges,
                error_weighted=weighted_errors,
                n_wavelength_bin=None,
            )
            if len(iq1d_elastic_wl) != 1:
                raise NotImplementedError("Not expected that there are more than 1 IQmod of elastic reference run.")

            iq2d_main_wl, iq1d_wl, k_vec, k_error_vec = normalize_by_elastic_reference_all(
                iq2d_main_wl,
                iq1d_main_wl[0],
                iq1d_elastic_wl[0],
                correction_setup.output_wavelength_dependent_profile,
                elastic_output_dir,
            )
            iq1d_main_wl[0] = iq1d_wl

            # write
            save_k_vector(
                iq1d_wl.wavelength,
                k_vec,
                k_error_vec,
                path=os.path.join(elastic_output_dir, f"{output_filename}_elastic_k1d_{k_file_prefix}.dat"),
            )

    # Inelastic incoherence correction
    if correction_setup.do_inelastic_correction[frameskip_frame]:
        inelastic_output_dir = os.path.join(
            output_dir, "info", "inelastic_incoh", output_filename, slice_name, f"frame_{frameskip_frame}"
        )
        os.makedirs(inelastic_output_dir, exist_ok=True)

        b_file_prefix = f"{raw_name}"

        # 1D correction
        iq2d_main_wl, iq1d_wl = do_inelastic_incoherence_correction(
            iq2d_main_wl,
            iq1d_main_wl[0],
            frameskip_frame,
            correction_setup,
            b_file_prefix,
            inelastic_output_dir,
            output_filename,
        )
        iq1d_main_wl[0] = iq1d_wl

    if not correction_setup.do_elastic_correction and not any(correction_setup.do_inelastic_correction):
        finite_iq1d = iq1d_in_frames[frameskip_frame]
        finite_iq2d = iq2d_in_frames[frameskip_frame]
        qmin = user_qmin
        qmax = user_qmax
    else:
        # Be finite
        finite_iq1d = iq1d_main_wl[0].be_finite()
        finite_iq2d = iq2d_main_wl.be_finite()
        # Bin binned I(Q1D, wl) and and binned I(Q2D, wl) in wavelength space
        assert len(iq1d_main_wl) == 1, (
            f"It is assumed that output I(Q) list contains 1 I(Q) but not {len(iq1d_main_wl)}"
        )

    # Bin output in wavelength space
    iq2d_main_out, iq1d_main_out = bin_all(
        finite_iq2d,
        finite_iq1d,
        num_x_bins,
        num_y_bins,
        n1dbins=num_q1d_bins,
        n1dbins_per_decade=num_q1d_bins_per_decade,
        decade_on_center=decade_on_center,
        bin1d_type=bin1d_type,
        log_scale=log_binning,
        qmin=qmin,
        qmax=qmax,
        qxrange=None,
        qyrange=None,
        annular_angle_bin=annular_bin,
        wedges=wedges,
        symmetric_wedges=symmetric_wedges,
        # When set to true, reduces high uncertainty in the high-Q limit when low statistics
        error_weighted=weighted_errors,
    )

    return iq2d_main_out, iq1d_main_out


def remove_workspaces(
    reduction_config: Dict,
    instrument_name: str,
    prefix: str,
    sample_run_number,
    center_run_number,
    extra_run_numbers: List,
):
    """Helping method to remove existing workspaces"""
    from drtsans.instruments import extract_run_number  # noqa E402
    from drtsans.path import registered_workspace  # noqa E402

    # In the future this should be made optional
    ws_to_remove = [f"{prefix}_{instrument_name}_{run_number}_raw_histo" for run_number in extra_run_numbers]
    # List special workspaces and workspace groups
    ws_to_remove.append(f"{prefix}_{instrument_name}_{sample_run_number}_raw_histo_slice_group")
    ws_to_remove.append(f"{prefix}_{instrument_name}_{center_run_number}_raw_events")
    ws_to_remove.append(f"{prefix}_sensitivity")
    ws_to_remove.append(f"{prefix}_mask")
    if "darkFileName" in reduction_config and reduction_config["darkFileName"]:
        run_number = extract_run_number(reduction_config["darkFileName"])
        ws_to_remove.append(f"{prefix}_{instrument_name}_{run_number}_raw_histo")
    if "blockedBeamRunNumber" in reduction_config and reduction_config["blockedBeamRunNumber"]:
        run_number = extract_run_number(reduction_config["blockedBeamRunNumber"])
        ws_to_remove.append(f"{prefix}_{instrument_name}_{run_number}_raw_histo")
    for ws_name in ws_to_remove:
        # Remove existing workspaces, this is to guarantee that all the data is loaded correctly
        if registered_workspace(ws_name):
            mtd.remove(ws_name)
