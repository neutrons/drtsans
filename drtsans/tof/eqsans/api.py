""" Top-level API for EQSANS """
from mantid.api import mtd
# Import rolled up to complete a single top-level API
from drtsans import (apply_sensitivity_correction, solid_angle_correction)
from drtsans import subtract_background
from drtsans.beam_finder import center_detector, find_beam_center
from drtsans.path import exists as path_exists
from drtsans.process_uncertainties import set_init_uncertainties  # noqa: F401
from drtsans.save_ascii import save_ascii_1D, save_xml_1D
from drtsans.save_2d import save_nist_dat, save_nexus
from drtsans.tof.eqsans.correct_frame import smash_monitor_spikes, transform_to_wavelength
from drtsans.tof.eqsans.load import load_events, load_events_monitor
from drtsans.tof.eqsans.dark_current import subtract_dark_current
from drtsans.mask_utils import apply_mask
from drtsans.tof.eqsans.normalization import normalize_by_flux

__all__ = ['apply_solid_angle_correction', 'subtract_background',
           'prepare_monitors', 'prepare_data', 'save_ascii_1D', 'save_xml_1D',
           'save_nist_dat', 'save_nexus', 'set_init_uncertainties']


def apply_solid_angle_correction(input_workspace):
    """Apply solid angle correction. This uses :func:`drtsans.solid_angle_correction`."""
    return solid_angle_correction(input_workspace,
                                  detector_type='VerticalTube')


def prepare_monitors(data, bin_width=0.1, output_workspace=None):
    r"""
    Loads monitor counts, correct TOF, and transforms to wavelength.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    bin_width: float
        Bin width for the output workspace, in Angstroms.
    output_workspace: str
        Name of the output workspace. If None, then it will be
        ``EQSANS_XXXXX_monitors`` with number XXXXX determined from ``data``.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    w = load_events_monitor(data, output_workspace=output_workspace)
    w = smash_monitor_spikes(w)
    w = transform_to_wavelength(w, bin_width=bin_width)
    return w


def prepare_data(data,
                 detector_offset=0, sample_offset=0,
                 bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                 center_x=None, center_y=None,
                 dark_current=None,
                 flux_method=None, flux=None,
                 mask=None, mask_panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None,
                 output_workspace=None, output_suffix=''):
    r"""
    Load an EQSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
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
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'proton charge',
        'monitor', or 'time'.
    flux: str
        if ``flux_method`` is proton charge, then path to file containing the
        wavelength distribution of the neutron flux. If ``flux method`` is
        monitor, then path to file containing the flux-to-monitor ratios.
        if ``flux_method`` is time, then pass one log entry name such
        as ``duration`` or leave it as :py:obj:`None` for automatic log search.
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
        Name of the output workspace. If not supplied, will be determined from the supplied value of ``data``.
    output_suffix: str
        If the ``output_workspace`` is not specified, this is appended to the automatically generated
        output workspace name.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    # First, load the event stream data into a workspace
    # The output_workspace name is for the Mantid workspace
    output_workspace = load_events(data, detector_offset=detector_offset,
                                   sample_offset=sample_offset,
                                   output_workspace=output_workspace, output_suffix=output_suffix)

    # The beam center should be provided by the reduction
    # script calling this function, but if it is not specified...
    if center_x is None or center_y is None:
        center_x, center_y = find_beam_center(output_workspace, mask=mask)
    center_detector(output_workspace, center_x=center_x, center_y=center_y)  # operates in-place

    # now that the instrument geometry has been set, we can convert
    # the event stream into neutrons binned into wavelength.
    # The conversion corrects for the different path lengths from
    # the sample to the detector pixels.
    # Interestingly, this function cannot handle constant dLambda/Lamba binning...
    output_workspace = transform_to_wavelength(output_workspace, bin_width=bin_width,
                                               low_tof_clip=low_tof_clip,
                                               high_tof_clip=high_tof_clip)

    # The initial estimate of the uncertainties are made next.
    output_workspace = set_init_uncertainties(output_workspace)

    # Next, we subtract dark current, if it exists.
    # Note that the function handles the normalization internally.
    if dark_current is not None:
        output_workspace = subtract_dark_current(output_workspace, dark_current)

    # The solid angle is corrected for next
    if solid_angle is True:
        output_workspace = apply_solid_angle_correction(output_workspace)

    # Interestingly, this is the only use of the btp dictionary.
    # The BTP stands for banks, tubes and pixels - it is a Mantid thing.
    apply_mask(output_workspace, panel=mask_panel, mask=mask, **btp)  # returns the mask

    # Correct for the detector sensitivity (the per pixel relative response)
    if sensitivity_file_path is not None \
            and path_exists(sensitivity_file_path):
        kw = dict(sensitivity_filename=sensitivity_file_path)
        output_workspace = apply_sensitivity_correction(output_workspace, **kw)

    # We can perform the desired normalization here.
    if flux_method is not None:
        kw = dict(method=flux_method)
        if flux_method == 'monitor':
            monitor_workspace = output_workspace + '_monitors'
            prepare_monitors(data, bin_width=bin_width,
                             output_workspace=monitor_workspace)
            kw['monitor_workspace='] = monitor_workspace
        output_workspace = normalize_by_flux(output_workspace, flux, **kw)

    if isinstance(output_workspace, str):
        return mtd[output_workspace]  # shouldn't happen
    else:
        return output_workspace
