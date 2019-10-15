from mantid.simpleapi import Subtract, mtd, DeleteWorkspace, Scale
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs

from drtsans.dark_current import duration

__all__ = ['subtract_dark_current', 'normalize_dark_current']


def normalize_dark_current(dark_workspace, output_workspace=None):
    r"""
    Divide a dark current workspace by its duration.

    Entry 'normalizing_duration' is added to the logs of the normalized
    dark current to annotate what log entry was used to find the duration

    **Mantid algorithms used:**
    :ref:`Scale <algm-Scale-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,

    Parameters
    ----------
    dark_workspace: str, ~mantid.api.MatrixWorkspace
        Dark current workspace
    output_workspace : str
        Name of the normalised dark workspace. If None, the name of the input
        workspace `dark_workspace` is chosen (and the input workspace is overwritten).

    Returns
    -------
    MatrixWorkspace
        Output workspace
    """
    if output_workspace is None:
        output_workspace = str(dark_workspace)

    # Find out the duration of the dark current from the logs, and divide
    dark_duration = duration(dark_workspace)
    normalizing_factor = unique_workspace_dundername()  # temporary workspace
    Scale(InputWorkspace=dark_workspace, Factor=1. / dark_duration.value, Operation='Multiply',
          OutputWorkspace=output_workspace)

    # Save the name of the log used to calculate the duration
    SampleLogs(output_workspace).insert('normalizing_duration', dark_duration.log_key)

    DeleteWorkspace(normalizing_factor)
    return mtd[output_workspace]


def subtract_dark_current(data_workspace, dark_workspace, output_workspace=None):
    r"""
    Subtract normalized dark from data, taking into account the duration of both the data and dark runs.

    ``normalized_data = data - (data_duration / dark_duration) * dark``

    **Mantid algorithms used:**
    :ref:`Scale <algm-Scale-v1>`,
    :ref:`Subtract <algm-Subtract-v1>`
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`

    Parameters
    ----------
    data_workspace: MatrixWorkspace
        Sample scattering with intensities versus wavelength.
    dark_workspace: MatrixWorkspace
        dark current workspace
    output_workspace : str
        Name of the output workspace. If None, the name of the input
        workspace `data_workspace` is chosen (and the input workspace is overwritten).

    Returns
    -------
    MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(data_workspace)

    # Normalize the dark current
    normalized_dark_current = unique_workspace_dundername()  # temporary workspace
    normalize_dark_current(dark_workspace, output_workspace=normalized_dark_current)

    # Find the duration of the data run using the same log key than that of the dark current
    duration_log_key = SampleLogs(normalized_dark_current).normalizing_duration.value
    data_duration = duration(data_workspace, log_key=duration_log_key).value
    Scale(InputWorkspace=normalized_dark_current, Factor=data_duration, Operation='Multiply',
          OutputWorkspace=normalized_dark_current)
    Subtract(LHSWorkspace=data_workspace, RHSWorkspace=normalized_dark_current, OutputWorkspace=output_workspace)

    DeleteWorkspace(normalized_dark_current)  # some clean-up
    return mtd[output_workspace]
