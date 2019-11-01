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


def normalize_by_time(input_workspace, output_workspace=None):
    """Normalise by time
    Used to normalise dark current

    **Mantid algorithms used:**
    :ref:`Divide <algm-Divide-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace
    output_workspace: str
        Optional name of the output workspace. Default is to replace the input workspace
    """
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    timer = SampleLogs(input_workspace).timer.value  # seconds
    timer_workspace = CreateSingleValuedWorkspace(timer, OutputWorkspace=unique_workspace_dundername())
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=timer_workspace, OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def normalize_by_monitor(input_workspace, output_workspace=None, factor_is=1.e08):
    """Normalise by the monitor value

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
