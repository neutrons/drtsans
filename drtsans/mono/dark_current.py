from mantid.simpleapi import Multiply, Subtract, mtd, CreateSingleValuedWorkspace, DeleteWorkspace, Scale
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs

from drtsans.dark_current import duration

__all__ = ['subtract_dark_current', 'normalise_to_workspace']


def normalise_to_workspace(dark_workspace, data_workspace, output_workspace=None):
    r"""
    Scale and Rebin in wavelength a `dark` current workspace with information
    from a `data` workspace.

    Rescale and rebin to the `data` workspace according to:
        frame_width_clipped / (frame_width * n_bins * duration) * I_dc(x, y)
    Entry 'normalizing_duration' is added to the logs of the normalized
    dark current to annotate what log entry was used to find the duration

    Parameters
    ----------
    dark_workspace: str, ~mantid.api.MatrixWorkspace
        Dark current workspace with units in time-of-flight
    data_workspace: str, ~mantid.api.MatrixWorkspace
        Sample scattering with intensities versus wavelength
    output_workspace : str
        Name of the normalised dark workspace. If None, the name of the input
        workspace `dark_workspace` is chosen (and the input workspace is overwritten).
    Returns
    -------
    MatrixWorkspace
        Output workspace, dark current rebinned to wavelength and rescaled
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
    Subtract normalized dark from data, taking into account the duration
    of both `data` and `dark` runs.
    Note that both of these input workspaces should have been normalised
    in advance.

    Parameters
    ----------
    data_workspace: MatrixWorkspace
        Sample scattering with intensities versus wavelength.
        NOT NORMALISED!
    dark_workspace: MatrixWorkspace
        Normalized dark current
    output_workspace : str
        Name of the output workspace. If None, the name of the input
        workspace `data_workspace` is chosen (and the input workspace is overwritten).
    Returns
    -------
    MatrixWorkspace
        `data` minus `dark` current
    """
    if output_workspace is None:
        output_workspace = str(data_workspace)

    # Normalize the dark current
    normalized_dark_current = unique_workspace_dundername()  # temporary workspace
    normalise_to_workspace(dark_workspace, data_workspace, output_workspace=normalized_dark_current)

    # Find the duration of the data run using the same log key than that of the dark current
    duration_log_key = SampleLogs(normalized_dark_current).normalizing_duration.value
    data_duration = duration(data_workspace, log_key=duration_log_key).value
    Scale(InputWorkspace=normalized_dark_current, Factor=data_duration, Operation='Multiply',
          OutputWorkspace=normalized_dark_current)
    Subtract(LHSWorkspace=data_workspace, RHSWorkspace=normalized_dark_current, OutputWorkspace=output_workspace)

    DeleteWorkspace(normalized_dark_current)  # some clean-up
    return mtd[output_workspace]
