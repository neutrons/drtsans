r"""
Links to mantid algorithms
CreateSingleValuedWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html>
Divide <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>
"""

from mantid.simpleapi import CreateSingleValuedWorkspace, Divide
from mantid.api import mtd

r"""
Hyperlinks to drtsans functions
SampleLogs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
"""  # noqa: E501
from drtsans.samplelogs import SampleLogs

__all__ = [
    "normalize_by_time",
    "normalize_by_monitor",
    "normalize_by_flux",
    "ZeroMonitorCountsError",
    "NoMonitorMetadataError",
]


class ZeroMonitorCountsError(ValueError):
    def __init__(self, message="Zero monitor counts found in the sample logs"):
        super().__init__(message)


class NoMonitorMetadataError(RuntimeError):
    def __init__(self, message="No monitor metadata found in the sample logs"):
        super().__init__(message)


def normalize_by_flux(ws, method):
    """Normalize to time or monitor

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    method : str
        Normalization method: 'time' or 'monitor'

    """
    if method == "time":
        return normalize_by_time(input_workspace=ws)
    elif method == "monitor":
        return normalize_by_monitor(input_workspace=ws)
    else:
        raise NotImplementedError()


def normalize_by_time(input_workspace, output_workspace=None):
    """Normalize by time
    Used to normalize dark current

    **Mantid algorithms used:**
    :ref:`Divide <algm-Divide-v1>`,

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    output_workspace: str
        Optional name of the output workspace. Default is to replace the input workspace
    """
    log_keys = ("duration", "timer")  # valid log keys to search for run duration
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    for log_key in log_keys:
        try:
            duration = SampleLogs(input_workspace).single_value(log_key)
            # Cast the timer value into a Mantid workspace to later divide the input workspace by this workspace
            duration_workspace = CreateSingleValuedWorkspace(duration, OutputWorkspace=mtd.unique_hidden_name())
            break
        except RuntimeError:
            continue  # check next log entry
    Divide(
        LHSWorkspace=input_workspace,
        RHSWorkspace=duration_workspace,
        OutputWorkspace=output_workspace,
    )
    duration_workspace.delete()  # some cleanup
    return mtd[output_workspace]


def normalize_by_monitor(input_workspace, output_workspace=None):
    """Normalize by the monitor value

    **Mantid algorithms used:**
    :ref: `CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>
    :ref:`Divide <algm-Divide-v1>`,

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    output_workspace: str
        Optional name of the output workspace. Default is to replace the input workspace

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        the normalized input workspace

    Raises
    ------
    RuntimeError
        No monitor metadata found in the sample logs of the input workspace
    """
    metadata_entry_names = [
        "monitor",  # created by load_events
        "monitor1",
    ]  # may be in the DAS logs

    reference_total_counts = 1.0e08  # actual number selected by the instrument team
    input_workspace = str(input_workspace)

    if output_workspace is None:
        output_workspace = input_workspace
    for entry_name in metadata_entry_names:
        monitor = None
        try:
            monitor = SampleLogs(input_workspace).single_value(entry_name) / reference_total_counts
            if float(monitor) <= 0.0:
                records = SampleLogs(input_workspace).records_json(logs=["run_start", "slice_info", "slice_interval"])
                raise ZeroMonitorCountsError(f"Zero monitor counts for workspace: {input_workspace}\n{records}")
            break
        except RuntimeError:  # the entry is not found in the metadata
            continue  # search next entry
    else:
        raise NoMonitorMetadataError(
            f"No monitor metadata found in the sample logs of the workspace: {input_workspace}"
        )
    # Cast the monitor value into a Mantid workspace to later divide the input workspace by this workspace
    monitor_workspace = CreateSingleValuedWorkspace(monitor, OutputWorkspace=mtd.unique_hidden_name())
    Divide(
        LHSWorkspace=input_workspace,
        RHSWorkspace=monitor_workspace,
        OutputWorkspace=output_workspace,
    )
    monitor_workspace.delete()
    return mtd[output_workspace]
