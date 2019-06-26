from __future__ import (absolute_import, division, print_function)
from ornl.sans.samplelogs import SampleLogsReader


def subtract_normalised_dark(data_ws, dark_current_ws):
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

    Returns
    -------
    MatrixWorkspace
        `data` minus `dark` current
    """

    duration_data = SampleLogsReader(data_ws).timer.value  # seconds
    __difference_normalised_dc = data_ws - duration_data * dark_current_ws
    return __difference_normalised_dc
