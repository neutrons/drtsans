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
        <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    output_workspace: str
        Optional name of the output workspace. Default is to replace the input workspace
    """
    log_keys = ('timer', 'duration')  # valid log keys to search for run duration
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    for log_key in log_keys:
        try:
            duration = SampleLogs(input_workspace).single_value('log_key')
            duration_workspace = CreateSingleValuedWorkspace(duration, OutputWorkspace=unique_workspace_dundername())
        except AttributeError:
            pass
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=duration_workspace, OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def normalize_by_monitor(input_workspace, output_workspace=None, factor_is=1.e08):
    """Normalize by the monitor value

    **Mantid algorithms used:**
    :ref: `CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>
         <https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html>
    :ref:`Divide <algm-Divide-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    output_workspace: str
        Optional name of the output workspace. Default is to replace the input workspace
    factor_is : float
        reference number of monitor counts selected by the instrument staff. Default is 10**8
    """
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    monitor = SampleLogs(input_workspace).monitor.value / factor_is  # seconds  # counts
    monitor_workspace = CreateSingleValuedWorkspace(monitor, OutputWorkspace=unique_workspace_dundername())
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=monitor_workspace, OutputWorkspace=output_workspace)
    return mtd[output_workspace]
