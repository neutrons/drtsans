from mantid.simpleapi import Multiply, Subtract, mtd, CreateSingleValuedWorkspace
from drtsans.settings import (unique_workspace_dundername as uwd)
from drtsans.dark_current import duration

__all__ = ['subtract_dark_current', 'normalise_to_workspace']


def normalise_to_workspace(dark_ws, data_ws, output_workspace=None):
    r"""
    Scale and Rebin in wavelength a `dark` current workspace with information
    from a `data` workspace.

    Rescale and rebin to the `data` workspace according to:
        frame_width_clipped / (frame_width * n_bins * duration) * I_dc(x, y)
    Entry 'normalizing_duration' is added to the logs of the normalized
    dark current to annotate what log entry was used to find the duration

    Parameters
    ----------
    dark_ws: str, ~mantid.api.MatrixWorkspace
        Dark current workspace with units in time-of-flight
    data_ws: str, ~mantid.api.MatrixWorkspace
        Sample scattering with intensities versus wavelength
    output_workspace : str
        Name of the normalised dark workspace. If None, the name of the input
        workspace `dark_ws` is chosen (and the input workspace is overwritten).
    Returns
    -------
    MatrixWorkspace
        Output workspace, dark current rebinned to wavelength and rescaled
    """
    if output_workspace is None:
        output_workspace = str(dark_ws)
    dark = mtd[str(dark_ws)]
    data = mtd[str(data_ws)]
    t_sam = duration(data).value
    t_data = duration(dark).value
    time_workspace = CreateSingleValuedWorkspace(t_sam/t_data, OutputWorkspace=uwd())
    Multiply(LHSWorkspace=time_workspace,
             RHSWorkspace=dark_ws,
             OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def subtract_dark_current(data_ws, dark_current_ws, output_workspace=None):
    r"""
    Subtract normalized dark from data, taking into account the duration
    of both `data` and `dark` runs.
    Note that both of these input workspaces should have been normalised
    in advance.

    Parameters
    ----------
    data_ws: MatrixWorkspace
        Sample scattering with intensities versus wavelength.
        NOT NORMALISED!
    dark_current_ws: MatrixWorkspace
        Normalized dark current
    output_workspace : str
        Name of the output workspace. If None, the name of the input
        workspace `data_ws` is chosen (and the input workspace is overwritten).
    Returns
    -------
    MatrixWorkspace
        `data` minus `dark` current
    """

    normalise_to_workspace(data_ws, dark_current_ws, output_workspace=str(dark_current_ws))
    Subtract(LHSWorkspace=data_ws,
             RHSWorkspace=dark_current_ws,
             OutputWorkspace=output_workspace)
    return mtd[output_workspace]
