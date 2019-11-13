""" BIOSANS API """
from mantid import simpleapi as msapi

from drtsans.mono import biosans
import drtsans
from drtsans.mask_utils import apply_mask
from drtsans.mono.load import load_events, transform_to_wavelength
from drtsans.mono.normalization import normalize_by_monitor, normalize_by_time
from drtsans.mono.dark_current import subtract_dark_current

# Functions exposed to the general user (public) API
__all__ = ['prepare_data']


def prepare_data(data,
                 mask_detector=None,
                 detector_offset=0, sample_offset=0,
                 center_x=None, center_y=None, center_y_wing=None,
                 dark_current=None,
                 flux_method=None,
                 mask=None, mask_panel=None, btp=dict(),
                 solid_angle=True,
                 sensitivity_file_path=None,
                 output_workspace=None, **kwargs):
    r"""
    Load a BIOSANS data file and bring the data to a point where it can be used. This includes applying basic
    corrections that are always applied regardless of whether the data is background or scattering data.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    mask_detector: str
        Name of an instrument component to mask
    detector_offset: float
        Additional translation of the detector along Z-axis, in mili-meters.
    sample_offset: float
        Additional translation of the sample along the Z-axis, in mili-meters.
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    center_y_wing: float
        Move the center of the wing detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
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
        Name of the output workspace.
    """
    # TODO: missing detector_offset and sample_offset
    ws = load_events(data, overwrite_instrument=True, output_workspace=output_workspace)
    ws = transform_to_wavelength(ws)

    if center_x is not None and center_y is not None and center_y_wing is not None:
        biosans.center_detector(ws, center_x=center_x,
                                center_y=center_y,
                                center_y_wing=center_y_wing)

    # Mask either detector
    if mask_detector is not None:
        msapi.MaskDetectors(ws, ComponentList=mask_detector)

    # Dark current
    if dark_current is not None:
        dark_ws = load_events(dark_current, overwrite_instrument=True)
        dark_ws = transform_to_wavelength(dark_ws)
        subtract_dark_current(ws, dark_ws)

    # Normalization
    if str(flux_method).lower() == 'monitor':
        normalize_by_monitor(ws)
    elif str(flux_method).lower() == 'time':
        normalize_by_time(ws)

    # Additional masks
    apply_mask(ws, panel=mask_panel, mask=mask, **btp)

    # Solid angle
    if solid_angle:
        biosans.solid_angle_correction(ws)

    # Sensitivity
    if sensitivity_file_path is not None:
        drtsans.apply_sensitivity_correction(ws, sensitivity_filename=sensitivity_file_path)

    return ws