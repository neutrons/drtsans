
def time(input_ws):
    """Normalise by time

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
    monitor = run.getProperty("monitor").value  # seconds

    res = input_ws / monitor
    return res
