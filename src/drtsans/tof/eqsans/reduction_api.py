# Move part of the methods from api.py to avoid importing in loops
import os
from collections import namedtuple
from typing import Dict, List, Tuple

import numpy as np
from mantid.simpleapi import (
    SaveAscii,
    mtd,
)  # noqa E402
from mantid.kernel import logger  # noqa E402

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
from drtsans.tof.eqsans.blocked_beam import subtract_blocked_beam  # noqa E402
from drtsans.tof.eqsans.correction_api import (
    CorrectionConfiguration,
)
from drtsans.tof.eqsans.dark_current import subtract_dark_current  # noqa E402
from drtsans.tof.eqsans.elastic_correction import (
    elastic_correction,
)
from drtsans.tof.eqsans.inelastic_correction import inelastic_correction
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
    blocked_beam=None,  # blocked beam workspace
    output_workspace=None,
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
    blocked_beam: ~mantid.dataobjects.Workspace2D
        histogram workspace containing the blocked beam measurement
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
        kw = dict(method=flux_method, output_workspace=output_workspace)
        if flux_method == "monitor":
            kw["monitor_workspace"] = data.monitor
        normalize_by_flux(output_workspace, flux, **kw)

    # Blocked Beam
    subtract_blocked_beam(
        output_workspace, blocked_beam, flux_method=flux_method, flux=flux, dark_current=dark_current
    )

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
        Objects containing 1D unbinned data I(|Q|). It will be used for scalar binned data
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
        Objects containing 1D unbinned data I(|Q|) for the elastic reference run
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
    # Get unbinned data for this frame
    iq2d = iq2d_in_frames[frameskip_frame]
    iq1d = iq1d_in_frames[frameskip_frame]

    # Preserve original unbinned data - both corrections should calculate factors from the SAME original data
    # This matches the old working pattern where one temporary binning was used for both corrections
    iq2d_for_factor_calc = iq2d
    iq1d_for_factor_calc = iq1d

    # Apply elastic correction if requested
    if correction_setup.do_elastic_correction and iq1d_elastic_ref_fr and iq2d_elastic_ref_fr:
        # Build output directory with slice and frame info
        elastic_dir = os.path.join(
            output_dir, "info", "elastic_norm", output_filename, slice_name, f"frame_{frameskip_frame}"
        )

        iq2d, iq1d = elastic_correction(
            iq2d_unbinned=iq2d,
            iq1d_unbinned=iq1d,
            iq2d_elastic_ref=iq2d_elastic_ref_fr[frameskip_frame],
            iq1d_elastic_ref=iq1d_elastic_ref_fr[frameskip_frame],
            num_x_bins=num_x_bins,
            num_y_bins=num_y_bins,
            num_q1d_bins=num_q1d_bins,
            num_q1d_bins_per_decade=num_q1d_bins_per_decade,
            decade_on_center=decade_on_center,
            bin1d_type=bin1d_type,
            log_binning=log_binning,
            user_qmin=user_qmin,
            user_qmax=user_qmax,
            annular_bin=annular_bin,
            wedges=wedges,
            symmetric_wedges=symmetric_wedges,
            weighted_errors=weighted_errors,
            output_wavelength_profile=correction_setup.output_wavelength_dependent_profile,
            output_dir=elastic_dir,
            output_filename=output_filename,
            raw_name=raw_name,
        )

    # Apply inelastic correction if requested
    if correction_setup.do_inelastic_correction[frameskip_frame]:
        # Build output directory with slice and frame info
        inelastic_dir = os.path.join(
            output_dir, "info", "inelastic_incoh", output_filename, slice_name, f"frame_{frameskip_frame}"
        )

        if correction_setup.do_elastic_correction and iq1d_elastic_ref_fr and iq2d_elastic_ref_fr:
            # Elastic correction was applied, so calculate b(λ) from original but apply to corrected
            from drtsans.iq import bin_all
            from drtsans.tof.eqsans.inelastic_correction import (
                calculate_incoherence_correction_factors,
                apply_incoherence_correction_to_unbinned_data,
            )

            # Determine Q ranges from original data (same as elastic correction would have used)
            qmin = user_qmin if user_qmin is not None else iq1d_for_factor_calc.mod_q.min()
            qmax = user_qmax if user_qmax is not None else iq1d_for_factor_calc.mod_q.max()
            qxrange = (np.min(iq2d_for_factor_calc.qx), np.max(iq2d_for_factor_calc.qx))
            qyrange = (np.min(iq2d_for_factor_calc.qy), np.max(iq2d_for_factor_calc.qy))

            # Temporarily bin ORIGINAL data for factor calculation
            _, iq1d_temp = bin_all(
                iq2d_for_factor_calc,
                iq1d_for_factor_calc,
                num_x_bins,
                num_y_bins,
                n1dbins=num_q1d_bins,
                n1dbins_per_decade=num_q1d_bins_per_decade,
                decade_on_center=decade_on_center,
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

            # Calculate b(λ) from temporarily binned ORIGINAL data
            correction_factors = calculate_incoherence_correction_factors(
                iq1d_temp[0],
                correction_setup.select_min_incoherence,
                correction_setup.select_intensityweighted[frameskip_frame],
                correction_setup.qmin[frameskip_frame],
                correction_setup.qmax[frameskip_frame],
                correction_setup.factor[frameskip_frame],
            )

            # Save b(λ)
            from drtsans.tof.eqsans.correction_api import save_b_factor
            from drtsans.tof.eqsans.inelastic_correction import CorrectedI1D
            from collections import namedtuple

            os.makedirs(inelastic_dir, exist_ok=True)
            WavelengthContainer = namedtuple("WavelengthContainer", ["wavelength"])
            wl_container = WavelengthContainer(wavelength=correction_factors.wavelength)
            save_b_factor(
                CorrectedI1D(wl_container, correction_factors.b_factor, correction_factors.b_error),
                os.path.join(inelastic_dir, f"{output_filename}_inelastic_b1d_{raw_name}.dat"),
            )

            # Apply b(λ) to ELASTIC-CORRECTED unbinned data
            logger.notice("Applying inelastic/incoherent correction to unbinned data")
            iq2d, iq1d = apply_incoherence_correction_to_unbinned_data(iq2d, iq1d, correction_factors)
        else:
            # No elastic correction, use the standard inelastic correction wrapper
            iq2d, iq1d = inelastic_correction(
                iq2d_unbinned=iq2d,
                iq1d_unbinned=iq1d,
                num_x_bins=num_x_bins,
                num_y_bins=num_y_bins,
                num_q1d_bins=num_q1d_bins,
                num_q1d_bins_per_decade=num_q1d_bins_per_decade,
                decade_on_center=decade_on_center,
                bin1d_type=bin1d_type,
                log_binning=log_binning,
                user_qmin=user_qmin,
                user_qmax=user_qmax,
                annular_bin=annular_bin,
                wedges=wedges,
                symmetric_wedges=symmetric_wedges,
                weighted_errors=weighted_errors,
                select_min_incoherence=correction_setup.select_min_incoherence,
                intensity_weighted=correction_setup.select_intensityweighted[frameskip_frame],
                incoh_qmin=correction_setup.qmin[frameskip_frame],
                incoh_qmax=correction_setup.qmax[frameskip_frame],
                incoh_factor=correction_setup.factor[frameskip_frame],
                output_dir=inelastic_dir,
                output_filename=output_filename,
                raw_name=raw_name,
            )

    # Remove non-finite values before final binning
    finite_iq2d = iq2d.be_finite()
    finite_iq1d = iq1d.be_finite()

    # ONE FINAL BINNING: Bin corrected (or uncorrected) unbinned data
    # This is the ONLY binning that produces output - the "One Rebin Only" pattern
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
        qmin=user_qmin,
        qmax=user_qmax,
        qxrange=None,
        qyrange=None,
        annular_angle_bin=annular_bin,
        wedges=wedges,
        symmetric_wedges=symmetric_wedges,
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
