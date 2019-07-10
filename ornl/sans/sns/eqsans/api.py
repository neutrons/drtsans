""" Top-level API for EQSANS """
# Import rolled up to complete a single top-level API
from .beam_finder import direct_beam_center
from ornl.settings import optional_output_workspace
from ornl.sans import solid_angle_correction

# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  center_detector, subtract_dark_current,
                                  normalise_by_flux, apply_mask)


def find_beam_center(ws, mask_file_path=None):
    """
        Beam center finder, with optional mask.
        input_ws: EventsWorkspace
            Workspace for the direct beam run
        mask_file_path: str
            File path of mask to apply before calculation

    """
    if mask_file_path is not None:
        ws = apply_mask(ws, mask_file_path)
    return direct_beam_center(ws)


def set_instrument_geometry(ws):
    """ Move detector components to their proper location """
    raise NotImplementedError()


def initial_uncertainty_estimation(ws):
    """
        Assign uncertainty for every bin in x, y, TOF
    """
    raise NotImplementedError()


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


def normalize(ws, normalization_type):
    """ Normalize to time, monitor, or proton charge """
    raise NotImplementedError()


@optional_output_workspace
def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=0.0, y_center=0.0,
                 dark_current=None,
                 flux=None,
                 mask=None, panel=None, btp=dict(),
                 sensitivity_file_path=None):
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
    ws = load_events(data, detector_offset=detector_offset,
                     sample_offset=sample_offset)
    ws = transform_to_wavelength(ws, bin_width=bin_width,
                                 low_tof_clip=low_tof_clip,
                                 high_tof_clip=high_tof_clip)
    center_detector(ws, x=x_center, y=y_center)
    if dark_current is not None:
        ws = subtract_dark_current(ws, dark_current)
    if flux is not None:
        ws = normalise_by_flux(ws, flux)
    apply_mask(ws, panel=panel, mask=mask, **btp)
    # Uncomment as we address them
    # ws = initial_uncertainty_estimation(ws)
    ws = apply_solid_angle_correction(ws)
    # ws = apply_sensitivity_correction(ws, sensitivity_file_path)
    return ws
