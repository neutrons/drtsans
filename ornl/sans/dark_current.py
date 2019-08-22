from __future__ import (absolute_import, division, print_function)

from dateutil.parser import parse as parse_date
import numpy as np
from mantid.simpleapi import (mtd, Integration, Transpose)
from ornl.settings import (namedtuplefy,
                           unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs


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
