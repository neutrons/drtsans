from drtsans.samplelogs import SampleLogs


def time(input_ws):
    """Normalise by time
    Used to normalise dark current

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """
    timer = SampleLogs(input_ws).timer.value  # seconds
    __time_normalisation = input_ws / timer
    return __time_normalisation


def monitor(input_ws):
    """Normalise by the monitor value

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """
    monitor = SampleLogs(input_ws).monitor.value  # seconds  # counts
    __monitor_normalisation = input_ws / monitor
    return __monitor_normalisation
