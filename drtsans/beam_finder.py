r""" Links to mantid algorithms
FindCenterOfMassPosition <https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v1.html>
Integration              <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
MoveInstrumentComponent  <https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html>
DeleteWorkspace          <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html>
"""
from mantid.simpleapi import FindCenterOfMassPosition, Integration, MoveInstrumentComponent, DeleteWorkspace
from mantid.kernel import logger

# drtsans imports
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.mask_utils import apply_mask, mask_spectra_with_special_values
from drtsans.solid_angle import solid_angle_correction

__all__ = ['center_detector', 'find_beam_center']  # exports to the drtsans namespace


def find_beam_center(input_workspace, method='center_of_mass', mask=None, mask_options={},
                     centering_options={}, solid_angle_method='VerticalTube'):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    **Mantid algorithms used:**
        :ref:`FindCenterOfMassPosition <algm-FindCenterOfMassPosition-v2>`,
        :ref:`Integration <algm-Integration-v1>`,

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
    method: str
        Method to calculate the beam center. Available methods are:
        - 'center_of_mass', invokes :ref:`FindCenterOfMassPosition <algm-FindCenterOfMassPosition-v1>`.
    mask: mask file path, `MaskWorkspace``, :py:obj:`list`.
        Mask to be passed on to ~drtsans.mask_utils.mask_apply.
    mask_options: dict
        Additional arguments to be passed on to ~drtsans.mask_utils.mask_apply.
    centering_options: dict
        Arguments to be passed on to the centering method.
    solid_angle_method: bool, str, specify which solid angle correction is needed

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters).
    """
    if method != 'center_of_mass':
        raise NotImplementedError()  # (f'{method} is not implemented')

    # integrate the TOF
    flat_ws = Integration(InputWorkspace=input_workspace, OutputWorkspace=uwd())
    mask_spectra_with_special_values(flat_ws)

    if mask is not None or mask_options != {}:
        mask_workspace = apply_mask(flat_ws, mask=mask, **mask_options)
        mask_workspace.delete()  # we don't need the mask workspace so keep it clean

    if solid_angle_method:
        solid_angle_correction(flat_ws, detector_type=solid_angle_method)

    # find center of mass position
    center = FindCenterOfMassPosition(InputWorkspace=flat_ws, **centering_options)
    logger.information("Found beam position: X={:.3} m, Y={:.3} m.".format(*center))
    DeleteWorkspace(flat_ws)
    return center


def center_detector(input_workspace, center_x, center_y, component='detector1'):
    """Translate the beam center currently located at (center_x, center_y) by an amount
    (-center_x, -center_y), so that the beam center is relocated to the origin of coordinates on the XY-plane

    **Mantid algorithms used:**
    :ref:`MoveInstrumentComponent <algm-MoveInstrumentComponent-v1>`,

    Parameters
    ----------
    input_workspace : Workspace2D, str
        The workspace to be centered
    center_x : float
        in meters
    center_y : float
        in meters
    component : string
        name of the detector to be centered
    """
    MoveInstrumentComponent(Workspace=input_workspace, ComponentName=component, X=-center_x, Y=-center_y,
                            RelativePosition=True)
