"""Helpers for reading EQSANS-specific sample logs."""

from typing import Union

import numpy as np
from mantid.api import MatrixWorkspace, Run

from drtsans.chopper import PULSE_FREQUENCY
from drtsans.samplelogs import SampleLogs


DEFAULT_FREQUENCY = PULSE_FREQUENCY
FREQUENCY_LOG = "frequency"
FALLBACK_FREQUENCY_LOG = "BL6:Det:TH:BL:Frequency"


def get_frequency(input_workspace: Union[str, Run, MatrixWorkspace, SampleLogs]) -> float:
    r"""Return the EQSANS pulse frequency in Hz.

    The primary ``frequency`` log is used when it contains a positive value.
    Some runs instead contain the detector timing log
    ``BL6:Det:TH:BL:Frequency``; its positive mean is used when the primary
    log is missing or zero.  A 60 Hz default is returned when neither log is
    usable.

    Parameters
    ----------
    input_workspace
        Workspace, run, workspace name, file accepted by :class:`SampleLogs`, or
        an existing :class:`SampleLogs` instance.

    Returns
    -------
    float
        Mean pulse frequency in Hz.
    """
    sample_logs = input_workspace if isinstance(input_workspace, SampleLogs) else SampleLogs(input_workspace)
    for log_name in (FREQUENCY_LOG, FALLBACK_FREQUENCY_LOG):
        try:
            frequency = sample_logs.single_value(log_name)
        except (AssertionError, AttributeError, KeyError, RuntimeError):
            continue
        if np.isfinite(frequency) and frequency > 0.0:
            return frequency
    return DEFAULT_FREQUENCY
