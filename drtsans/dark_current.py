from dateutil.parser import parse as parse_date
import numpy as np

# Links to mantid algorithms
# Integration <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
# DeleteWorkspace <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html>
from mantid.simpleapi import mtd, Integration, DeleteWorkspace

from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.samplelogs import SampleLogs


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
    Find the total number of neutron counts in each detector pixel.

    In a detector pixel has no counts, then the error of the zero counts is set to one.

    Parameters
    ----------
    input_workspace: str, EventsWorkspace
        Usually a dark current workspace for which we need to know the total number of counts per pixel-detector

    Returns
    -------
    tuple
        counts, error in the counts
    """
    # Create a workspace containing the total counts per pixel, and starting errors
    counts_workspace = unique_workspace_dundername()
    Integration(input_workspace, OutputWorkspace=counts_workspace)

    counts = mtd[counts_workspace].extractY().flatten()
    errors = mtd[counts_workspace].extractE().flatten()
    errors[np.where(counts == 0)[0]] = 1

    DeleteWorkspace(counts_workspace)  # some clean-up
    return counts, errors
