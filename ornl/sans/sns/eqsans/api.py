""" Top-level API for EQSANS """
from mantid.api import mtd
from mantid import simpleapi as sapi
from ornl.settings import unique_workspace_dundername as uwd
# Import rolled up to complete a single top-level API
from ornl.sans import apply_sensitivity_correction, solid_angle_correction
# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  center_detector, subtract_dark_current,
                                  normalise_by_flux, apply_mask)
from ornl.path import exists as path_exists

__all__ = ['apply_solid_angle_correction', 'subtract_background',
           'prepare_data']


def apply_solid_angle_correction(ws):
    """
        Apply solid angle correction
    """
    return solid_angle_correction(ws, detector_type='VerticalTube')


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


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=None, y_center=None,
                 dark_current=None,
                 flux=None,
                 mask=None, panel=None, btp=dict(),
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
    flux: str
        Path to file containing the wavelength distribution
        of the neutron flux.
    panel: str
        Either 'front' or 'back' to mask a whole panel
    mask: mask file path, MaskWorkspace
        Mask to be applied.
    btp: dict
        Additional properties to MaskBTP, if desired
    """
    # let load_events dictate the name of the workspace
    output_workspace = load_events(data, detector_offset=detector_offset,
                                   sample_offset=sample_offset,
                                   output_workspace=output_workspace)
    output_workspace = str(output_workspace)  # convert it to its name
    transform_to_wavelength(output_workspace, bin_width=bin_width,
                            low_tof_clip=low_tof_clip,
                            high_tof_clip=high_tof_clip)
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)
    if flux is not None:
        normalise_by_flux(output_workspace, flux)
    apply_mask(output_workspace, panel=panel, mask=mask, **btp)
    apply_solid_angle_correction(output_workspace)
    center_detector(output_workspace, x=x_center, y=y_center)
    if sensitivity_file_path is not None \
            and path_exists(sensitivity_file_path):
        apply_sensitivity_correction(output_workspace,
                                     sensitivity_filename=sensitivity_file_path)  # noqa: E501 can't fit in line
    return mtd[output_workspace]
