# This module contains legacy code to do sensitivity correction with detector patching method.
# It is not supposed to be used or supported after the new algorithm passes acceptance.
r"""
Links to mantid algorithms
https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html
https://docs.mantidproject.org/nightly/algorithms/SaveNexusProcessed-v1.html

https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
https://docs.mantidproject.org/nightly/algorithms/LoadNexusProcessed-v2.html
https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html
https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html
https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
"""
from mantid.simpleapi import mtd, CalculateEfficiency, MaskDetectorsIf, SaveNexusProcessed, \
    DeleteWorkspace, Divide, LoadNexusProcessed, MaskDetectors, \
    Integration, CreateWorkspace
from mantid.kernel import Property, logger

__all__ = ['calculate_sensitivity_correction']


# This is old sensitivity calculation method.  It is kept temporarily for comparison purpose
def calculate_sensitivity_correction(input_workspace, min_threshold=0.5, max_threshold=2.0,
                                     filename=None, output_workspace=None):
    """
    Calculate the detector sensitivity

    **Mantid algorithms used:**
    :ref:`CalculateEfficiency <algm-CalculateEfficiency-v1>`,
    :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
    :ref:`SaveNexusProcessed <algm-SaveNexusProcessed-v1>`


    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace to calculate the sensitivity from
    min_threshold: float
        Minimum threshold for efficiency value.
    max_threshold: float
        Maximum threshold for efficiency value
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    """
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    CalculateEfficiency(InputWorkspace=input_workspace, OutputWorkspace=output_workspace,
                        MinThreshold=min_threshold, MaxThreshold=max_threshold)
    MaskDetectorsIf(InputWorkspace=output_workspace, OutputWorkspace=output_workspace,
                    Mode='SelectIf', Operator='Equal', Value=Property.EMPTY_DBL)
    if filename is not None:
        SaveNexusProcessed(InputWorkspace=output_workspace, Filename=filename)
    return mtd[output_workspace]
