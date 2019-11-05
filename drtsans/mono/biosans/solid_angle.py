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

    # The first two spectra could contain monitor detectors, e.g. if using LoadEmptyInstrument or if using
    # LoadInstrument(RewriteSpectraMap=True)
    shift = 2 if mtd[str(input_workspace)].spectrumInfo().isMonitor(0) is True else 0

    # Apply correction to main detector
    spectra_in_main_detector = 24 * 8 * 256  # there are 24 eight-packs, 8 tubes/eight-pack, 256 pixels/tube
    first_index, last_index = shift, shift + spectra_in_main_detector - 1
    solid_angle.solid_angle_correction(input_workspace, detector_type='VerticalTube',
                                       StartWorkspaceIndex=first_index, EndWorkspaceIndex=last_index,
                                       output_workspace=output_workspace)
    # Apply correction to wing detector
    spectra_in_wing_detector = 20 * 8 * 256  # there are 20 eight-packs
    first_index = shift + spectra_in_main_detector
    last_index = shift + spectra_in_main_detector + spectra_in_wing_detector - 1
    solid_angle.solid_angle_correction(output_workspace, detector_type='VerticalWing',
                                       StartWorkspaceIndex=first_index, EndWorkspaceIndex=last_index)
    return mtd[output_workspace]
