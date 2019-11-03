r"""
Links to Mantid algorithms
DeleteWorkspace https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
Divide https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
SolidAngle https://docs.mantidproject.org/nightly/algorithms/SolidAngle-v1.html
ReplaceSpecialValues https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html
"""
from mantid.api import mtd

r"""
Links to drtsans modules and functions
solid_angle <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/solid_angle.py>
"""  # noqa: E501
from drtsans import solid_angle

__all__ = ['solid_angle_correction', ]


def solid_angle_correction(input_workspace, output_workspace=None):
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
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of the input workspace is taken,
        thus the output workspace replaces the input workspace.

    Returns
    -------
    ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    solid_angle.solid_angle_correction(input_workspace, detector_type='VerticalTube', StartWorkspaceIndex=2,
                                       EndWorkspaceIndex=49153, output_workspace=output_workspace)
    # Apply correction to wing detector
    solid_angle.solid_angle_correction(output_workspace, detector_type='VerticalWing', StartWorkspaceIndex=49154,
                                       EndWorkspaceIndex=90113)
    return mtd[output_workspace]
