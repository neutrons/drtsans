from __future__ import print_function


def time(input_ws):
    """Normalise by time
    Used to normalise dark current

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """

    run = input_ws.getRun()
    timer = run.getProperty("timer").value  # seconds

    res = input_ws / timer
    return res


def monitor(input_ws):
    """Normalise by the monitor value

    Parameters
    ----------
    input_ws : [Mantid Workspace]

    """

    run = input_ws.getRun()
    monitor = run.getProperty("monitor").value  # counts

    res = input_ws / monitor
    return res
