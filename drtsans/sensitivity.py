import os
import numpy as np
from mantid.kernel import Property, logger
from drtsans.settings import unique_workspace_name as uwn
from drtsans.path import exists as path_exists
r"""
Links to mantid algorithms
https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/CalculateEfficiency-v1.html
https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
https://docs.mantidproject.org/nightly/algorithms/LoadNexusProcessed-v2.html
https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html
https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html
https://docs.mantidproject.org/nightly/algorithms/SaveNexusProcessed-v1.html
https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html
https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
"""
from mantid.simpleapi import mtd, CloneWorkspace, CalculateEfficiency, \
    DeleteWorkspace, Divide, LoadNexusProcessed, MaskDetectors, \
    MaskDetectorsIf, SaveNexusProcessed, \
    Integration, CreateWorkspace

__all__ = ['apply_sensitivity_correction', 'calculate_sensitivity_correction']


# flake8: noqa: C901
def apply_sensitivity_correction(input_workspace, sensitivity_filename=None,
                                 sensitivity_workspace=None,
                                 min_threshold=None,
                                 max_threshold=None,
                                 output_workspace=None):
    """Apply a previously calculated sensitivity correction

    **Mantid algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,
    :ref:`Divide <algm-Divide-v1>`,
    :ref:`LoadNexusProcessed <algm-LoadNexusProcessed-v1>`,
    :ref:`MaskDetectors <algm-MaskDetectors-v1>`
    :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        workspace to apply the correction to
    sensitivity_filename: str
        file containing previously calculated sensitivity correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    min_threshold: float or None
        if not None, the data will be masked if the sensitivity
        is below this threshold
    max_threshold: float or None
        if not None, the data will be masked if the sensitivity
        is above this threshold
    output_workspace:  ~mantid.api.MatrixWorkspace
        corrected workspace. This is the input workspace by default
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    cleanupSensitivity = False
    if sensitivity_workspace is None or \
            str(sensitivity_workspace) not in mtd:  # load the file
        if sensitivity_workspace is None:
            sensitivity_workspace = os.path.split(sensitivity_filename)[-1]
            sensitivity_workspace = sensitivity_workspace.split('.')[0]
        cleanupSensitivity = True
        if not path_exists(sensitivity_filename):
            msg = 'Cannot find file "{}"'.format(sensitivity_filename)
            raise RuntimeError(msg)
        LoadNexusProcessed(Filename=sensitivity_filename,
                           OutputWorkspace=sensitivity_workspace)
    if (not sensitivity_workspace) or (str(sensitivity_workspace) not in mtd):
        raise RuntimeError('No sensitivity workspace provided')

    if str(input_workspace) != str(output_workspace):
        CloneWorkspace(InputWorkspace=input_workspace,
                       OutputWorkspace=output_workspace)
    MaskDetectors(Workspace=output_workspace,
                  MaskedWorkspace=sensitivity_workspace)

    # additional masking dependent on threshold
    temp_sensitivity = CloneWorkspace(InputWorkspace=sensitivity_workspace,
                                      OutputWorkspace=uwn(prefix="__sensitivity_"))

    # Process masked pixels in the sensitivities workspace, i.e., pixels with value NaN
    min_threshold = process_masked_pixels(temp_sensitivity, min_threshold)

    if min_threshold is not None:
        MaskDetectorsIf(InputWorkspace=temp_sensitivity,
                        Operator='LessEqual',
                        Value=min_threshold,
                        OutputWorkspace=temp_sensitivity)
    if max_threshold is not None:
        MaskDetectorsIf(InputWorkspace=temp_sensitivity,
                        Operator='GreaterEqual',
                        Value=max_threshold,
                        OutputWorkspace=temp_sensitivity)
    Divide(LHSWorkspace=output_workspace, RHSWorkspace=temp_sensitivity,
           OutputWorkspace=output_workspace)
    DeleteWorkspace(temp_sensitivity)

    if cleanupSensitivity:
        DeleteWorkspace(sensitivity_workspace)

    # set empty units
    mtd[output_workspace].setYUnit('')

    return mtd[output_workspace]


def process_masked_pixels(sensitivity_workspace, min_threshold):
    """Convert data value of masked pixels (NaN) to a value less than min_threshold, for example, 0

    If min_threshold is None, then set the  min_threshold to a reasonable small value, i.e., 1E-6 and

    Parameters
    ----------
    sensitivity_workspace :  ~mantid.api.MatrixWorkspace
        sensitivity workspace
    min_threshold : float or None
        minimum threshold

    Returns
    -------
    float
        minimum threshold (new)

    """
    # Process threshold
    if min_threshold is None:
        min_threshold = 1E-6
        nan_value = 0
    else:
        nan_value = min(0., min_threshold - 1E-6)

    # Loop through workspace
    num_spec = sensitivity_workspace.getNumberHistograms()
    for i_ws in range(num_spec):
        if np.isnan(sensitivity_workspace.readX(i_ws)[0]):
            sensitivity_workspace.dataX(i_ws)[0] = 1.
            sensitivity_workspace.dataY(i_ws)[0] = nan_value

    return min_threshold


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
