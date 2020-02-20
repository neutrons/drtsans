from mantid.simpleapi import ClearMaskFlag, DeleteWorkspace, Divide, mtd, SolidAngle, ReplaceSpecialValues
from drtsans.instruments import instrument_enum_name
from drtsans.settings import unique_workspace_dundername

__all__ = ['calculate_solid_angle', 'solid_angle_correction']


def calculate_solid_angle(input_workspace, detector_type='VerticalTube', output_workspace=None, **kwargs):
    '''Calculate the solid angle from the ``input_workspace``.

    **Mantid algorithms used:**
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
    '''
    instrument = instrument_enum_name(input_workspace)

    # default behaviour is to create instrument unique name
    if not output_workspace:
        output_workspace = '{}_solid_angle'.format(instrument)
    # return existing workspace if possible
    if output_workspace in mtd:
        return mtd[output_workspace]

    # calculate the solid angle
    SolidAngle(InputWorkspace=input_workspace, OutputWorkspace=output_workspace, Method=detector_type)

    if kwargs:  # assume the pixel range was set
        # set the solid angle of the mystery parts of the instrument as 1 and don't mask them
        ClearMaskFlag(Workspace=output_workspace)
        ReplaceSpecialValues(InputWorkspace=output_workspace, OutputWorkspace=output_workspace,
                             SmallNumberThreshold=1.e-9, SmallNumberValue=1.)

    return mtd[output_workspace]


def solid_angle_correction(input_workspace, detector_type='VerticalTube', output_workspace=None,
                           solid_angle_ws=None, **kwargs):
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
    solid_angle_ws: str, ~mantid.api.MatrixWorkspace
        Workspace containing the solid angle calculation
    kwargs: dict
        Additional arguments to Mantid algorithm :ref:`SolidAngle <algm-SolidAngle-v1>`

    Returns
    -------
    ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    """
    input_workspace = str(input_workspace)

    if output_workspace is None:
        output_workspace = input_workspace

    # calculate the solid angle for the instrument if the workspace wasn't provided
    cleanup_solid_angle_ws = False
    if not solid_angle_ws:
        cleanup_solid_angle_ws = True
        solid_angle_ws = calculate_solid_angle(input_workspace=input_workspace, detector_type=detector_type,
                                               output_workspace=unique_workspace_dundername(), **kwargs)

    # correct the input workspace and get rid of nan and infinity
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=solid_angle_ws, OutputWorkspace=output_workspace)
    ReplaceSpecialValues(InputWorkspace=output_workspace, NaNValue=0., InfinityValue=0.,
                         OutputWorkspace=output_workspace)

    # delete temporary workspaces
    if cleanup_solid_angle_ws:
        DeleteWorkspace(solid_angle_ws)

    return mtd[output_workspace]
