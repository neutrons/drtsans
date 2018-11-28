from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.kernel import logger
from mantid.simpleapi import (CloneWorkspace, Minus, Integration, Transpose,
                              RebinToWorkspace)


from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import bank_detector_ids
from ornl.sans.sns.eqsans.correct_frame import frame_width, tof_clippings


def compute_log_ratio(data, dark, log_key):
    """
    Compute ratio of data to dark duration for one log entry

    Parametes
    ---------
    data: MatrixWorkspace
        Scattering data workspace
    dark: MatrixWorkspace
        Dark current workspace
    log_key: str
        Name of the log containing the duration

    Returns
    -------
    float

    Raises
    ------
    NoneType
        Log entry is not present
    """
    dt = SampleLogs(data)[log_key]
    dk = SampleLogs(dark)[log_key]
    try:
        dt = dt.getStatistics().duration
        dk = dk.getStatistics().duration
    except AttributeError:
        dt = dt.value
        dk = dk.value
    return dt / dk


def duration_ratio(data, dark, log_key=None):
    """
    Compute the ratio of data to dark-current duration. Return 1.0 if the
    duration ration cannot be computed

    Parameters
    ----------
    data: MatrixWorkspace
        Scattering data workspace
    dark: MatrixWorkspace
        Dark current workspace
    log_key: None or str
        Name of the log containing the duration. If `None`, duration will
        be tried looking sequentially into log entries 'duration',
        'proton_charge', and 'timer

    Returns
    -------
    float
    """
    not_found = 'Logs could not be found, duration ratio set to 1.0'
    if log_key is not None:
        try:
            return compute_log_ratio(data, dark, log_key)
        except AttributeError:
            logger.error(not_found)
    else:
        for l in ('duration', 'proton_charge', 'timer'):
            try:
                return compute_log_ratio(data, dark, l)
            except AttributeError:
                continue
    return 1.0


def duration(dark, log_key=None):
    """
    Compute the dark-current duration. Return 1.0 if the
    duration ratio cannot be computed

    Parameters
    ----------
    dark: MatrixWorkspace
        Dark current workspace
    log_key: None or str
        Name of the log containing the duration. If `None`, duration will
        be tried looking sequentially into log entries 'duration',
        'proton_charge', and 'timer

    Returns
    -------
    float
    """
    pass


def counts_in_detector(dark):
    r"""
    Create a workspace with the total number of neutron counts in each detector

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
    _wnc = Integration(dark)
    _wnc = Transpose(_wnc)
    y = _wnc.dataY(0)  # counts
    e = _wnc.dataE(0)  # errors
    e[y == 0] = 1.0  # convention of error=1 if no counts present

    _wnc.setE(0, e)
    _wnc = Transpose(_wnc)
    return _wnc, y


def distribute_counts(dark, data):
    r"""
    Distribute the neutron counts over the wavelength channels of `data`

    Parameters
    ----------
    dark: Workspace2D
        Single bin histograms, containing the total neutron count per detector
    data: Workspace2D
        Data workspace, rebinned in wavelength
    """
    _dark = RebinToWorkspace(dark, data, OuptutWorkspace=dark.name())


def normalise_to_workspace(dark, data, out_ws, log_key=None):
    r"""
    Scale and rebin in wavelength a `dark` current workspace with information
    from a `data` workspace.

    First, normalize by the duration of the dark current run. Also
    rescale and rebin to the `data` workspace according to:
    $\frac{1}{dark_duration} \frac{frame_width-(low_tof_cut + high_tof_cut)}{frame_width} \frac{\Delta \lambda_j}{\lambda_{max} - \lamdba_{min}} I_{dc}(x, y)$

    Parameters
    ----------
    dark: EventsWorkspace
        Dark current workspace
    data: MatrixWorkspace
        Sample scattering
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
    """  # noqa
    _darkn, nc = counts_in_detector(dark)
    distribute_counts(_darkn, data)

        RebinToWorkspace(cid.ws, data)  # conform to wavelength binning

    ws = CloneWorkspace(data, OutputWorkspace=out_ws)
    x = ws.dataX(0)
    dx = x[1:] - x[:-1]  # \Delta \lambda_j
    fr = frame_width(data)
    ltc, htc = tof_clippings(data)
    dr = duration(dark, log_key=log_key)
    cid = counts_in_detector(dark)  # total neutron counts in each detector

    y = (1.0 / dr) * ((fr - (ltc + htc)) / fr) * dx / (x[-1] - x[0])
    e = y / np.sqrt(n)
    for i in bank_detector_ids(ws, masked=False):
        ws.dataY(i)[:] = y
    return ws


@namedtuplefy
def subtract_dark_current(data, dark, out_ws, out_darkn='darkn',
                           log_key=None):
    r"""

    Parameters
    ----------
    dark: EventsWorkspace
        Dark current workspace
    data: MatrixWorkspace
        Sample scattering
    out_ws: str
        Name of the output workspace containing scattering minus dark current
    out_darkn: str
        Name of the output workspace containing the normalized dark current
    log_key: str
        Log key to search for duration of the runs. if `None`, the function
        does a sequential search for entries 'duration', 'proton_charge',
        and 'timer'

    Returns
    -------
    nametuple
        Fields of the name tuple:
        - data: workspace containing scattering minus dark current
        - dark: workspace containing the normalized dark current
    """
    darkn = normalise_to_workspace(dark, data, out_ws=out_darkn, log_key=None)
    ws = Minus(data, darkn, OutputWorkspace=out_ws)
    return dict(data=ws, dark=darkn)