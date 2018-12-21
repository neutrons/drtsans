from __future__ import print_function
from ornl.sans.samplelogs import SampleLogs


def time(input_ws):
    """Normalise by time
    Used to normalise dark current

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """
    timer = SampleLogs(input_ws).timer.value  # seconds
    res = input_ws / timer
    return res


def monitor(input_ws):
    """Normalise by the monitor value

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """
    monitor = SampleLogs(input_ws).monitor.value  # seconds  # counts
    res = input_ws / monitor
    return res
