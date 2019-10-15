import numpy as np
from mantid.simpleapi import (mtd, RebinToWorkspace, ConvertUnits,
                              Subtract, Scale, LoadEventNexus)

from drtsans.dark_current import counts_in_detector, duration
from drtsans.settings import (amend_config,
                              unique_workspace_dundername as uwd)
from drtsans.path import exists, registered_workspace
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import correct_frame

__all__ = ['subtract_dark_current', 'normalize_dark_current']


def normalize_dark_current(dark_workspace, data_workspace, output_workspace=None):
    r"""
    Scale and Rebin in wavelength a ``dark`` current workspace with information
    from a ``data`` workspace.

    Rescale and rebin to the ``data`` workspace according to:

    .. math:: frame\_width\_clipped / (frame\_width * n\_bins * duration) * I\_dc(x, y)

    Entry 'normalizing_duration' is added to the logs of the normalized
    dark current to annotate what log entry was used to find the duration

    **Mantid algorithms used:**
    :ref:`ConvertUnits <algm-ConvertUnits-v1>`,
    :ref:`RebinToWorkspace <algm-RebinToWorkspace-v1>`

    Parameters
    ----------
    dark_workspace: str, EventsWorkspace
        Dark current workspace with units in time-of-flight
    data_workspace: str, MatrixWorkspace
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
    sample_logs = SampleLogs(data_workspace)
    # rescale counts by the shorter considered TOF width
    tof_clipping_factor = sample_logs.tof_frame_width_clipped.value / sample_logs.tof_frame_width.value
    counts_in_pixel, counts_error_in_pixel = counts_in_detector(dark_workspace)  # counts and error per detector
    ConvertUnits(InputWorkspace=dark_workspace, Target='Wavelength', Emode='Elastic', OutputWorkspace=output_workspace)
    RebinToWorkspace(WorkspaceToRebin=output_workspace, WorkspaceToMatch=data_workspace, PreserveEvents=False,
                     OutputWorkspace=output_workspace)
    #
    # Determine the histogram bins for which data should have counts
    # If running in frame-skipped mode, there is a wavelength band
    # gap between the lead and skipped bands
    #
    dark_normalized = mtd[output_workspace]
    bands = correct_frame.clipped_bands_from_logs(data_workspace)  # lead and pulse bands
    gap_indexes = correct_frame.band_gap_indexes(data_workspace, bands)
    n_gap_bins = len(gap_indexes)
    n_bins = len(dark_normalized.dataY(0))
    n_significant_bins = n_bins - n_gap_bins  # wavelength bins with counts
    #
    # factor_y is number of counts per unit time and wavelength bin
    #
    dark_duration = duration(dark_workspace)
    counts_rescaling_factor = tof_clipping_factor / (dark_duration.value * n_significant_bins)
    counts_error_rescaling_factor = tof_clipping_factor / (dark_duration.value * np.sqrt(n_significant_bins))
    #
    # Rescale the dark counts. Pay attention to the wavelength gap in frame-skip mode
    #
    for i in range(dark_normalized.getNumberHistograms()):
        dark_normalized.dataY(i)[:] = counts_rescaling_factor * counts_in_pixel[i]
        dark_normalized.dataE(i)[:] = counts_error_rescaling_factor * counts_error_in_pixel[i]
    if n_gap_bins > 0:
        for i in range(dark_normalized.getNumberHistograms()):
            dark_normalized.dataY(i)[gap_indexes] = 0.0
            dark_normalized.dataE(i)[gap_indexes] = 0.0
    SampleLogs(dark_normalized).insert('normalizing_duration', dark_duration.log_key)
    return dark_normalized


def subtract_normalised_dark_current(input_workspace, dark_ws,
                                     output_workspace=None):
    r"""
    Subtract normalized dark current from data, taking into account
    the duration of both 'data' and 'dark' runs.

    Entry 'normalizing_duration' is added to the logs of the output workspace
    to annotate what log entry was used to find the duration of both
    'data' and 'dark' runs. Log entry 'normalizing_duration' must be
    present in the logs of workspace 'dark'.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Sample scattering with intensities versus wavelength
    dark_ws: str, ~mantid.api.MatrixWorkspace
        Normalized dark current after being normalized with
        `normalize_dark_current`
    output_workspace : str
        Name of the workspace after dark current subtraction. If :py:obj:`None`,
        the name of the input workspace is chosen (and the input workspace
        is overwritten).

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        'data' minus 'dark' current
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    duration_log_key = SampleLogs(dark_ws).normalizing_duration.value
    d = duration(input_workspace, log_key=duration_log_key).value
    scaled = Scale(InputWorkspace=dark_ws, Factor=d, OutputWorkspace=uwd())
    Subtract(LHSWorkspace=input_workspace, RHSWorkspace=scaled, OutputWorkspace=output_workspace)
    scaled.delete()
    SampleLogs(output_workspace).insert('normalizing_duration', duration_log_key)
    return mtd[output_workspace]


def subtract_dark_current(input_workspace, dark, output_workspace=None):
    r"""


    Parameters
    ----------
    input_workspace : int, str, ~mantid.api.IEventWorkspace
        The workspace to be normalised
    dark: int, str, ~mantid.api.IEventWorkspace
        run number, file path, workspace name, or :py:obj:`~mantid.api.IEventWorkspace`
        for dark current.
    output_workspace : str
        Name of the workspace after dark current subtraction. If None,
        the name of the input workspace is chosen (and the input workspace
        is overwritten).

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    if registered_workspace(dark):
        _dark = dark
    elif (isinstance(dark, str) and exists(dark)) or isinstance(dark, int):
        with amend_config({'default.instrument': 'EQSANS'}):
            _dark = LoadEventNexus(Filename=dark, OutputWorkspace=uwd())
    else:
        message = 'Unable to find or load the dark current {}'.format(dark)
        raise RuntimeError(message)

    _dark_normal = normalize_dark_current(_dark, input_workspace,
                                          output_workspace=uwd())
    subtract_normalised_dark_current(input_workspace, _dark_normal,
                                     output_workspace=output_workspace)
    _dark_normal.delete()
    if _dark is not dark:
        _dark.delete()

    return mtd[output_workspace]
