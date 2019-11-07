r"""
Links to mantid algorithms
CreateSingleValuedWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html>
Divide <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>
"""
from mantid.simpleapi import CreateSingleValuedWorkspace, Divide
from mantid.api import mtd

r"""
Hyperlinks to drtsans functions
unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
SampleLogs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
"""  # noqa: E501
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs

__all__ = ['normalize_by_time', 'normalize_by_monitor']


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
    log_keys = ('duration', 'timer')  # valid log keys to search for run duration
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    for log_key in log_keys:
        try:
            duration = SampleLogs(input_workspace).single_value(log_key)
            # Cast the timer value into a Mantid workspace to later divide the input workspace by this workspace
            duration_workspace = CreateSingleValuedWorkspace(duration, OutputWorkspace=unique_workspace_dundername())
            break
        except RuntimeError:
            continue  # check next log entry
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=duration_workspace, OutputWorkspace=output_workspace)
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
    """
    reference_total_counts = 1.e08  # actual number selected by the instrument team
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    monitor = SampleLogs(input_workspace).monitor.value / reference_total_counts
    # Cast the monitor value into a Mantid workspace to later divide the input workspace by this workspace
    monitor_workspace = CreateSingleValuedWorkspace(monitor, OutputWorkspace=unique_workspace_dundername())
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=monitor_workspace, OutputWorkspace=output_workspace)
    return mtd[output_workspace]
