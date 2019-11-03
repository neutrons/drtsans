r"""
Links to Mantid algorithms
DeleteWorkspace https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
Divide https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
SolidAngle https://docs.mantidproject.org/nightly/algorithms/SolidAngle-v1.html
ReplaceSpecialValues https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html
"""
from mantid.simpleapi import DeleteWorkspace, Divide, mtd, SolidAngle, ReplaceSpecialValues
r"""
Links to drtsans functions
unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
"""  # noqa: E501
from drtsans.settings import unique_workspace_dundername

__all__ = ['solid_angle_correction']


def solid_angle_correction(input_workspace, detector_type='VerticalTube', output_workspace=None, **kwargs):
    r"""
    The algorithm calculates solid angles subtended by the individual pixel-detectors when vieved from the sample
    position. The returned workspace is the input workspace normalized (divided) by the pixel solid angles.

    **Mantid algorithms used:**
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
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=solid_angle_ws, OutputWorkspace=output_workspace)
    DeleteWorkspace(solid_angle_ws)
    ReplaceSpecialValues(InputWorkspace=output_workspace, NaNValue=0., InfinityValue=0.,
                         OutputWorkspace=output_workspace)
    return mtd[output_workspace]
