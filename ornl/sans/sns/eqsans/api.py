""" Top-level API for EQSANS """
from mantid.api import mtd
# Import rolled up to complete a single top-level API
from ornl.sans import solid_angle_correction
# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  center_detector, subtract_dark_current,
                                  normalise_by_flux, apply_mask)

__all__ = ['prepare_data']


def apply_solid_angle_correction(ws):
    """
        Apply solid angle correction
    """
    return solid_angle_correction(ws, detector_type='VerticalTube')


def apply_sensitivity_correction(ws, sensitivity_file_path):
    """
        Apply sensitivity correction
    """
    raise NotImplementedError()


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=0.0, y_center=0.0,
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
    center_detector(output_workspace, x=x_center, y=y_center)
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)
    if flux is not None:
        normalise_by_flux(output_workspace, flux)
    apply_mask(output_workspace, panel=panel, mask=mask, **btp)
    # Uncomment as we address them
    # initial_uncertainty_estimation(output_workspace)
    apply_solid_angle_correction(output_workspace)
    # apply_sensitivity_correction(output_workspace, sensitivity_file_path)
    return mtd[output_workspace]
