r""" Links to mantid algorithms
FindCenterOfMassPosition <https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v1.html>
Integration              <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
MoveInstrumentComponent  <https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html>
"""
from mantid.simpleapi import FindCenterOfMassPosition, Integration, MoveInstrumentComponent, CreateWorkspace, Divide
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
        pixel_areas = bounding_box_widths[:, 0] * bounding_box_widths[:, 2]
        number_Of_spectra = flat_ws.getNumberHistograms()
        X_axis_values = flat_ws.readX(0)
        workspace_pixelarea = CreateWorkspace(DataX=X_axis_values, DataY=pixel_areas,
                                              Nspec=number_Of_spectra, OutputWorkspace='area')

        area_corrected_counts = Divide(LHSWorkspace=flat_ws, RHSWorkspace=workspace_pixelarea,
                                       OutputWorkspace='corrected_counts')
        flat_ws = Integration(area_corrected_counts, OutputWorkspace=uwd())

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
