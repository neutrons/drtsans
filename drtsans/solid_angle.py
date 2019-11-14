from mantid.simpleapi import ClearMaskFlag, DeleteWorkspace, Divide, mtd, SolidAngle, ReplaceSpecialValues
from drtsans.settings import unique_workspace_dundername

__all__ = ['solid_angle_correction']


def solid_angle_correction(input_workspace, detector_type='VerticalTube', output_workspace=None, **kwargs):
    r"""
    The algorithm calculates solid angles subtended by the individual pixel-detectors when vieved from the sample
    position. The returned workspace is the input workspace normalized (divided) by the pixel solid angles.

    **Mantid algorithms used:**
    :ref:`ClearMaskFlag <algm-ClearMaskFlag-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,
    :ref:`Divide <algm-Divide-v1>`,
    :ref:`ReplaceSpecialValues <algm-ReplaceSpecialValues-v1>`,
    :ref:`SolidAngle <algm-SolidAngle-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Input workspace to be normalized by the solid angle.
    detector_type: str
        Select the method to calculate the Solid Angle. Allowed values: [‘GenericShape’,
        ‘Rectangle’, ‘VerticalTube’, ‘HorizontalTube’, ‘VerticalWing’, ‘HorizontalWing’]
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of the input workspace is taken,
        thus the output workspace replaces the input workspace.
    kwargs: dict
        Additional arguments to Mantid algorithm :ref:`SolidAngle <algm-SolidAngle-v1>`

    Returns
    -------
    ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    """
    input_workspace = str(input_workspace)
    solid_angle_ws = unique_workspace_dundername()  # temporary workspace
    if output_workspace is None:
        output_workspace = input_workspace

    SolidAngle(InputWorkspace=input_workspace, OutputWorkspace=solid_angle_ws, Method=detector_type, **kwargs)
    if kwargs:  # assume the pixel range was set
        # set the solid angle of the mystery parts of the instrument as 1 and don't mask them
        ClearMaskFlag(Workspace=solid_angle_ws)
        ReplaceSpecialValues(InputWorkspace=solid_angle_ws, OutputWorkspace=solid_angle_ws,
                             SmallNumberThreshold=1.e-9, SmallNumberValue=1.)

    # correct the input workspace and get rid of nan and infinity
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=solid_angle_ws, OutputWorkspace=output_workspace)
    DeleteWorkspace(solid_angle_ws)
    ReplaceSpecialValues(InputWorkspace=output_workspace, NaNValue=0., InfinityValue=0.,
                         OutputWorkspace=output_workspace)
    return mtd[output_workspace]
