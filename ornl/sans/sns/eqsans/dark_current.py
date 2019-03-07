from __future__ import (absolute_import, division, print_function)

from dateutil.parser import parse as parse_date
import numpy as np
from mantid.simpleapi import (Integration, Transpose, RebinToWorkspace,
                              ConvertUnits, RenameWorkspace)

from ornl.settings import unique_workspace_name, namedtuplefy
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans import correct_frame as cf


@namedtuplefy
def duration(dark, log_key=None):
    """
    Compute the duration of the workspace by iteratively searching the logs for
    keys 'duration', 'start_time/end_time', 'proton_charge', and 'timer'.

    Parameters
    ----------
    dark: MatrixWorkspace
        Workspace, usually the dark current workspace
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
    log_keys = ('duration', 'start_time', 'proton_charge', 'timer') if \
        log_key is None else (log_key, )
    sl = SampleLogs(dark)

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


def counts_in_detector(dark):
    r"""
    Fin the total number of neutron counts in each detector pixel.
    By definition, error=1 when zero counts in the detector.

    Parameters
    ----------
    dark: EventsWorkspace
        Deark current workspace

    Returns
    -------
    tuple
        Elements of the tuple
        - ws: Workspace2D, Single bin histograms, containing the total neutron
            count per detector
        - n: numpy.ndarrary, list of counts
    """
    _wnc = Integration(dark, OutputWorkspace=unique_workspace_name())
    _wnc = Transpose(_wnc, OutputWorkspace=_wnc.name())
    y = np.copy(_wnc.dataY(0))  # counts
    e = np.copy(_wnc.dataE(0))  # errors
    _wnc.delete()
    e[y < 1] = 1.0  # convention of error=1 if no counts present
    return y, e


def normalise_to_workspace(dark, data, out_ws):
    r"""
    Scale and rebin in wavelength a `dark` current workspace with information
    from a `data` workspace.

    Rescale and rebin to the `data` workspace according to:
        frame_width_clipped / (frame_width * n_bins * duration) * I_dc(x, y)
    Entry 'normalizing_duration' is added to the logs of the output workspace
    to annotate what log entry was used to find the duration

    Parameters
    ----------
    dark: EventsWorkspace
        Dark current workspace with units in time-of-flight
    data: MatrixWorkspace
        Sample scattering with intensities versus wavelength
    out_ws: str
        Name of the normalized output workspace
    log_key: str
        Log key to search for duration of the runs. if `None`, the function
        does a sequential search for entries 'duration', 'proton_charge',
        and 'timer'

    Returns
    -------
    MatrixWorkspace
        Output workspace, dark current rebinned to wavelength and rescaled
    """
    sl = SampleLogs(data)
    fwr = sl.tof_frame_width_clipped.value / sl.tof_frame_width.value
    nc, ec = counts_in_detector(dark)  # counts and error per detector
    _dark = ConvertUnits(dark, Target='Wavelength', Emode='Elastic',
                         OutputWorkspace=out_ws)
    _dark = RebinToWorkspace(_dark, data, PreserveEvents=False)
    bands = cf.clipped_bands_from_logs(data)  # lead and pulse bands
    gap_indexes = cf.band_gap_indexes(data, bands)
    n_gap_bins = len(gap_indexes)
    n_bins = len(_dark.dataY(0))
    n_significant_bins = n_bins - n_gap_bins
    d = duration(dark)
    factor = fwr / (d.value * n_significant_bins)
    for i in range(_dark.getNumberHistograms()):
        _dark.dataY(i)[:] = np.ones(n_bins) * factor * nc[i]
        _dark.dataE(i)[:] = np.ones(n_bins) * factor * ec[i]
    if n_gap_bins > 0:
        for i in range(_dark.getNumberHistograms()):
            _dark.dataY(i)[gap_indexes] = 0.0
            _dark.dataE(i)[gap_indexes] = 0.0
    SampleLogs(_dark).normalizing_duration = d.log_key  # append to the logs
    return _dark


def subtract_normalised_dark(data, dark, out_ws):
    r"""
    Subtract normalized dark from data, taking into account the duration
    of both `data` and `dark` runs.

    Entry 'normalizing_duration' is added to the logs of the output workspace
    to annotate what log entry was used to find the duration of both
    `data` and `dark` runs.

    Parameters
    ----------
    data: MatrixWorkspace
        Sample scattering with intensities versus wavelength
    dark: MatrixWorkspace
        Normalized dark current after being normalized with
        `normalise_to_workspace`
    out_ws: str
        Name of the subtracted output workspace

    Returns
    -------
    MatrixWorkspace
        `data` minus `dark` current
    """
    duration_log_key = SampleLogs(dark).normalizing_duration.value
    d = duration(data, log_key=duration_log_key).value
    difference = data - d * dark
    RenameWorkspace(difference, out_ws)
    SampleLogs(difference).normalizing_duration = duration_log_key
    return difference
