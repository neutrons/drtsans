""" Top-level API for EQSANS """
# Import rolled up to complete a single top-level API
from .beam_finder import direct_beam_center
from mantid.simpleapi import mtd
from ornl.sans import solid_angle_correction

# Imports from EQSANS public API
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  center_detector, subtract_dark_current,
                                  normalise_by_flux)


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
    return solid_angle_correction(ws, detector_type='VerticalTube')


def apply_sensitivity_correction(ws, sensitivity_file_path):
    """
        Apply sensitivity correction
    """
    raise NotImplementedError()


def normalize(ws, normalization_type):
    """ Normalize to time, monitor, or proton charge """
    raise NotImplementedError()


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 x_center=0.0, y_center=0.0,
                 dark_current=None,
                 flux=None,
                 mask_file_path=None,
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
    """
    # let load_events dictate the name of the workspace
    output_workspace = load_events(data, detector_offset=detector_offset,
                                   sample_offset=sample_offset,
                                   output_workspace=output_workspace)
    output_workspace = str(output_workspace)  # convert it to its name
    transform_to_wavelength(output_workspace, bin_width=bin_width,
                            low_tof_clip=low_tof_clip,
                            high_tof_clip=high_tof_clip,
                            output_workspace=output_workspace)
    center_detector(output_workspace, x=x_center, y=y_center)
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current,
                              output_workspace=output_workspace)
    if flux is not None:
        normalise_by_flux(output_workspace, flux,
                          output_workspace=output_workspace)
    # Uncomment as we address them
    # apply_mask(ws, mask_file_path, output_workspace=output_workspace)
    # initial_uncertainty_estimation(ws, output_workspace=output_workspace)
    apply_solid_angle_correction(output_workspace,
                                 output_workspace=output_workspace)
    # apply_sensitivity_correction(ws, sensitivity_file_path,
    #                              output_workspace=output_workspace)
    return mtd[output_workspace]
