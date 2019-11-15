r""" Links to mantid algorithms
FindCenterOfMassPosition <https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v1.html>
Integration              <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
MoveInstrumentComponent  <https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html>
CreateWorkspace          <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
DeleteWorkspace          <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html>
"""
from mantid.simpleapi import FindCenterOfMassPosition, Integration, MoveInstrumentComponent, CreateWorkspace, \
    DeleteWorkspace
from mantid.kernel import logger
import numpy as np

# drtsans imports
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.mask_utils import apply_mask, mask_spectra_with_special_values


__all__ = ['center_detector', 'find_beam_center']  # exports to the drtsans namespace


def find_beam_center(input_workspace, method='center_of_mass', mask=None, mask_options={},
                     centering_options={}, area_corection_flag=True):
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
    area_corection_flag: str, flag to specify if area correction is needed

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters).
    """
    if method != 'center_of_mass':
        raise NotImplementedError()  # (f'{method} is not implemented')

    # flatten the workspace
    flat_ws = Integration(InputWorkspace=input_workspace, OutputWorkspace=uwd())
    mask_spectra_with_special_values(flat_ws)

    if mask is not None or mask_options != {}:
        mask_workspace = apply_mask(flat_ws, mask=mask, **mask_options)
        mask_workspace.delete()  # we don't need the mask workspace so keep it clean

    if area_corection_flag:
        bounding_box_widths = np.array(
            [flat_ws.getDetector(i).shape().getBoundingBox().width() for i in
             range(flat_ws.getNumberHistograms())])
        # bounding_box_widths is an array of shape=(N, 3) where N is the number of detectors.
        # Thus, bounding_box_widths[:, 0] are the widths of the detectors along the first axis (X-axis),
        # bounding_box_widths[:, 1] are the widths of the detectors along the second axis (Y-axis) and
        # bounding_box_widths[:, 2] are the widths of the detectors along the third axis (Z-axis).
        # The pixel area should be the width along the X-axis multiplied by the width along the Y-axis.
        pixel_areas = bounding_box_widths[:, 0] * bounding_box_widths[:, 1]
        number_Of_spectra = flat_ws.getNumberHistograms()
        X_axis_values = flat_ws.readX(0)
        workspace_pixelarea = CreateWorkspace(DataX=X_axis_values, DataY=pixel_areas,
                                              Nspec=number_Of_spectra, OutputWorkspace='area')

        # Mantid allows to work with handles to workspaces or with the strings containing the names of the workspaces.
        # If using handles, the '/' operator is used. If using strings, Divide algorithm is used.
        # Here we need the handle, so we simply used '/' operator.
        flat_ws /= workspace_pixelarea
        DeleteWorkspace(workspace_pixelarea)

    # find center of mass position
    center = FindCenterOfMassPosition(InputWorkspace=flat_ws, **centering_options)
    logger.information("Found beam position: X={:.3} m, Y={:.3} m.".format(*center))
    return center


def center_detector(input_workspace, center_x, center_y, component='detector1'):
    """Center the detector. Move the `component` by (-center_x, -centery)
    from the current position (relative motion).

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
