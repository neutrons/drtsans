""" Top-level API for EQSANS """
# Import rolled up to complete a single top-level API
from .beam_finder import direct_beam_center
from ornl.settings import optional_output_workspace

# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  center_detector)


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


def divide_by_flux(ws, flux_file_path):
    """ Divide wavelength distribution by the flux distribution """
    raise NotImplementedError()


def apply_mask(ws, mask_file_path):
    """
        Apply mask to detector pixels
    """
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
    raise NotImplementedError()


def apply_sensitivity_correction(ws, sensitivity_file_path):
    """
        Apply sensitivity correction
    """
    raise NotImplementedError()


def normalize(ws, normalization_type):
    """ Normalize to time, monitor, or proton charge """
    raise NotImplementedError()


@optional_output_workspace
def prepare_data(file_path, detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=0.0, y_center=0.0,
                 dark_current=None,
                 mask_file_path=None, sensitivity_file_path=None):
    """
        Load an EQSANS data file and bring the data to a point where it
        can be used. This includes applying basic corrections that are
        always applied regardless of whether the data is background or
        scattering data.
    """
    ws = load_events(file_path, detector_offset=detector_offset,
                     sample_offset=sample_offset)
    ws = transform_to_wavelength(ws, bin_width=bin_width,
                                 low_tof_clip=low_tof_clip,
                                 high_tof_clip=high_tof_clip)
    center_detector(ws, x=x_center, y=y_center)
    # Uncomment as we address them

    # ws = subtract_dark_current(ws, dark_current)
    # ws = apply_mask(ws, mask_file_path)
    # ws = initial_uncertainty_estimation(ws)
    # ws = apply_solid_angle_correction(ws)
    # ws = apply_sensitivity_correction(ws, sensitivity_file_path)
    # ws = divide_by_flux(ws)
    # ws = normalize(ws)
    return ws