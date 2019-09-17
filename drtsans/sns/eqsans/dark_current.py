from __future__ import (absolute_import, division, print_function)

from dateutil.parser import parse as parse_date
import numpy as np
from mantid.simpleapi import (mtd, Integration, Transpose, RebinToWorkspace,
                              ConvertUnits, Subtract, Scale, LoadEventNexus)

from drtsans.settings import (namedtuplefy, amend_config,
                           unique_workspace_dundername as uwd)
from drtsans.path import exists, registered_workspace
from drtsans.samplelogs import SampleLogs
from drtsans.sns.eqsans import correct_frame

__all__ = ['subtract_dark_current', ]


@namedtuplefy
def duration(input_workspace, log_key=None):
    """
    Compute the duration of the workspace by iteratively searching the logs for
    keys 'duration', 'start_time/end_time', 'proton_charge', and 'timer'.

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace
        Usually the dark current workspace
    log_key: str
        If a log entry is passed, only the contents under this log entry are
        searched. No iterative search over the default values is performed.

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - value: float, contents under the log
        - log_key: str, log used to return the duration

    """
    ws = mtd[str(input_workspace)]
    log_keys = ('duration', 'start_time', 'proton_charge', 'timer') if \
        log_key is None else (log_key, )
    sl = SampleLogs(ws)

    def from_start_time(lk):
        st = parse_date(sl[lk].value)
        et = parse_date(sl['end_time'].value)
        return (et - st).total_seconds()

    def from_proton_charge(lk):
        return sl[lk].getStatistics().duration

    calc = dict(start_time=from_start_time, proton_charge=from_proton_charge)

    for lk in log_keys:
        try:
            return dict(value=calc.get(lk, sl.single_value)(lk), log_key=lk)
        except RuntimeError:
            continue
    raise AttributeError("Could not determine the duration of the run")


def counts_in_detector(input_workspace):
    r"""
    Fin the total number of neutron counts in each detector pixel.
    By definition, error=1 when zero counts in the detector.

    Parameters
    ----------
    input_workspace: str, EventsWorkspace
        Usually the dark current workspace

    Returns
    -------
    tuple
        Two elements in the tuple: (1)numpy.ndarray: counts;
        (2) numpy.ndarray: error in the counts
    """
    ws = mtd[str(input_workspace)]
    _wnc = Integration(ws, OutputWorkspace=uwd())
    _wnc = Transpose(_wnc, OutputWorkspace=_wnc.name())
    y = np.copy(_wnc.dataY(0))  # counts
    e = np.copy(_wnc.dataE(0))  # errors
    _wnc.delete()
    e[y < 1] = 1.0  # convention of error=1 if no counts present
    return y, e


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
    dark_ws: str, EventsWorkspace
        Dark current workspace with units in time-of-flight
    data_ws: str, MatrixWorkspace
        Sample scattering with intensities versus wavelength
    log_key: str
        Log key to search for duration of the runs. if `None`, the function
        does a sequential search for entries 'duration', 'proton_charge',
        and 'timer'
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
    sl = SampleLogs(data)
    # rescale counts by the shorter considered TOF width
    fwr = sl.tof_frame_width_clipped.value / sl.tof_frame_width.value
    nc, ec = counts_in_detector(dark)  # counts and error per detector
    ConvertUnits(InputWorkspace=dark,
                 Target='Wavelength', Emode='Elastic',
                 OutputWorkspace=output_workspace)
    RebinToWorkspace(WorkspaceToRebin=output_workspace, WorkspaceToMatch=data,
                     PreserveEvents=False,
                     OutputWorkspace=output_workspace)
    _dark = mtd[output_workspace]
    #
    # Determine the histogram bins for which data should have counts
    # If running in frame-skipped mode, there is a wavelength band
    # gap between the lead and skipped bands
    #
    bands = correct_frame.clipped_bands_from_logs(data)  # lead and pulse bands
    gap_indexes = correct_frame.band_gap_indexes(data, bands)
    n_gap_bins = len(gap_indexes)
    n_bins = len(_dark.dataY(0))
    n_significant_bins = n_bins - n_gap_bins  # wavelength bins with counts
    #
    # factor_y is number of counts per unit time and wavelength bin
    #
    dark_duration = duration(dark)
    factor_y = fwr / (dark_duration.value * n_significant_bins)
    factor_e = fwr / (dark_duration.value * np.sqrt(n_significant_bins))
    #
    #
    #
    for i in range(_dark.getNumberHistograms()):
        _dark.dataY(i)[:] = factor_y * nc[i]
        _dark.dataE(i)[:] = factor_e * ec[i]
    if n_gap_bins > 0:
        for i in range(_dark.getNumberHistograms()):
            _dark.dataY(i)[gap_indexes] = 0.0
            _dark.dataE(i)[gap_indexes] = 0.0
    SampleLogs(_dark).insert('normalizing_duration', dark_duration.log_key)
    return _dark


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
        `normalise_to_workspace`
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
    scaled = Scale(InputWorkspace=dark_ws, Factor=-d, OutputWorkspace=uwd())
    Subtract(LHSWorkspace=input_workspace, RHSWorkspace=scaled,
             OutputWorkspace=output_workspace)
    scaled.delete()
    SampleLogs(output_workspace).insert('normalizing_duration',
                                        duration_log_key)
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

    _dark_normal = normalise_to_workspace(_dark, input_workspace,
                                          output_workspace=uwd())
    subtract_normalised_dark_current(input_workspace, _dark_normal,
                                     output_workspace=output_workspace)
    _dark_normal.delete()
    if _dark is not dark:
        _dark.delete()

    return mtd[output_workspace]
