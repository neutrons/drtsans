# https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
from mantid.simpleapi import CreateSingleValuedWorkspace, Divide, mtd
from drtsans.samplelogs import SampleLogs
from drtsans.settings import (unique_workspace_dundername as uwd)


def time(input_workspace, output_workspace=None):
    """Normalise by time
    Used to normalise dark current

    Parameters
    ----------
    input_workspace : [Mantid Workspace]

    """
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace
    timer = SampleLogs(input_workspace).timer.value  # seconds
    timer_workspace = CreateSingleValuedWorkspace(timer, OutputWorkspace=uwd())
    Divide(LHSWorkspace=input_workspace,
           RHSWorkspace=timer_workspace,
           OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def monitor(input_workspace, output_workspace=None, factor_is=10 ** 8):
    """Normalise by the monitor value

    **Mantid algorithms used:**
    :ref: `CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>
    :ref:`Divide <algm-Divide-v1>`,

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
    monitor_workspace = CreateSingleValuedWorkspace(monitor, OutputWorkspace=uwd())
    Divide(LHSWorkspace=input_workspace,
           RHSWorkspace=monitor_workspace,
           OutputWorkspace=output_workspace)
    return mtd[output_workspace]
