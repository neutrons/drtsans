""" Top-level API for EQSANS """
from mantid.api import mtd
from mantid import simpleapi as sapi
from ornl.settings import unique_workspace_dundername as uwd
# Import rolled up to complete a single top-level API
from ornl.sans import (apply_sensitivity_correction, solid_angle_correction)
from ornl.sans.save_ascii import save_ascii_1D, save_xml_1D
from ornl.sans.save_2d import save_nist_dat, save_nexus
from ornl.sans.process_uncertainties import set_init_uncertainties  # noqa: F401
# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, load_events_monitor,
                                  transform_to_wavelength,
                                  center_detector, subtract_dark_current,
                                  normalise_by_flux, apply_mask)
from ornl.sans.sns.eqsans.correct_frame import smash_monitor_spikes
from ornl.sans.thickness_normalization import normalize_by_thickness
from ornl.path import exists as path_exists

__all__ = ['apply_solid_angle_correction', 'subtract_background',
           'prepare_monitors', 'prepare_data', 'save_ascii_1D', 'save_xml_1D',
           'save_nist_dat', 'save_nexus']


def apply_solid_angle_correction(input_workspace):
    """
        Apply solid angle correction
    """
    return solid_angle_correction(input_workspace,
                                  detector_type='VerticalTube')


def normalize(ws, normalization_type):
    """ Normalize to time, monitor, or proton charge """
    raise NotImplementedError()


def subtract_background(input_workspace, background, scale=1.0,
                        output_wokspace=None):
    r"""
    Subtract a prepared background from a prepared sample.

    Perform a rebin if sample and background have different binning.

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace
        Sample workspace.
    background: str, MatrixWorkspace
        Background workspace.
    scale: float
        Rescale background intensities by this multiplicative factor before
        subtraction from the sample.
    output_wokspace: str
        Name of the sample corrected by the background. If None, then
        `input_workspace` will be overwritten.

    Returns
    -------
    MatrixWorkspace
    """
    if output_wokspace is None:
        output_wokspace = str(input_workspace)
    ws = mtd[str(input_workspace)]
    wb = mtd[str(background)]
    wb2 = sapi.RebinToWorkspace(WorkspaceToRebin=wb,
                                WorkspaceToMatch=ws,
                                OutputWorkspace=uwd())
    wb2 = sapi.Scale(InputWorkspace=wb2, OutputWorkspace=wb2.name(),
                     Factor=scale, Operation='Multiply')
    sapi.Minus(LHSWorkspace=ws, RHSWorkspace=wb2,
               OutputWorkspace=output_wokspace)
    wb2.delete()
    return mtd[output_wokspace]


def prepare_monitors(data, bin_width=0.1, output_workspace=None):
    r"""
    Loads monitor counts, correct TOF, and transforms to wavelength.

    Parameters
    ----------
    data: int, str, EventWorkspace
        Run number as int or str, file path, EventWorkspace
    bin_width: float
        Bin width for the output workspace, in Angstroms.
    output_workspace: str
        Name of the output workspace. If None, then it will be
        EQSANS_XXXXX_monitors with number XXXXX determined from `data`.

    Returns
    -------
    MatrixWorkspace
    """
    w = load_events_monitor(data, output_workspace=output_workspace)
    w = smash_monitor_spikes(w)
    w = transform_to_wavelength(w, bin_width=bin_width)
    return w


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=None, y_center=None,
                 dark_current=None,
                 flux_method=None, flux=None,
                 mask=None, mask_panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None,
                 output_workspace=None):
    r"""
    Load an EQSANS data file and bring the data to a point where it
    can be used. This includes applying basic corrections that are
    always applied regardless of whether the data is background or
    scattering data.

    Parameters
    ----------
    data: int, str, EventWorkspace
        Run number as int or str, file path, EventWorkspace
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
    x_center: float
        Move the center of the detector to this X-coordinate. If `None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have `x=0`.
    y_center: float
        Move the center of the detector to this X-coordinate. If `None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have `x=0`.
    dark_current: int, str, EventWorkspace
        Run number as int or str, file path, EventWorkspace
    flux_method: str
        Method for flux normalization. Either 'proton charge',
        'monitor', or 'time'.
    flux: str
        if `flux_method` is proton charge, then path to file containing the
        wavelength distribution of the neutron flux. If `flux method` is
        monitor, then path to file containing the flux-to-monitor ratios.
        if `flux_method` is time, then pass one log entry name such
        as `duration` or leave it as `None` for automatic log search.
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
    output_workspace: str
        Name of the output workspace. If None, then it will be
        EQSANS_XXXXX with number XXXXX determined from the supplied `data`.
    """
    # let load_events dictate the name of the workspace
    output_workspace = load_events(data, detector_offset=detector_offset,
                                   sample_offset=sample_offset,
                                   output_workspace=output_workspace)
    output_workspace = str(output_workspace)  # convert it to its name
    transform_to_wavelength(output_workspace, bin_width=bin_width,
                            low_tof_clip=low_tof_clip,
                            high_tof_clip=high_tof_clip)
    center_detector(output_workspace, x=x_center, y=y_center)
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)
    # Normalization by flux
    if flux_method is not None:
        kw = dict(method=flux_method)
        if flux_method == 'monitor':
            monitor_workspace = output_workspace + '_monitors'
            prepare_monitors(data, bin_width=bin_width,
                             output_workspace=monitor_workspace)
            kw['monitor_workspace='] = monitor_workspace
        normalise_by_flux(output_workspace, flux, **kw)
    apply_mask(output_workspace, panel=mask_panel, mask=mask, **btp)
    if solid_angle is True:
        apply_solid_angle_correction(output_workspace)
    if sensitivity_file_path is not None \
            and path_exists(sensitivity_file_path):
        kw = dict(sensitivity_filename=sensitivity_file_path)
        apply_sensitivity_correction(output_workspace, **kw)
    return mtd[output_workspace]
