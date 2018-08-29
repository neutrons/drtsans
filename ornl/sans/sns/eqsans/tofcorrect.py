from __future__ import (absolute_import, division, print_function)


def frame_width(workspace, pulse_freq='frequency'):
    """
    Time between pulses, in micro-seconds

    Parameters
    ----------
    workspace:  Mantid.Workspace

    pulse_freq: str
        Entry name in the log containing the pulses generation frequencies

    Returns
    -------
    float

    Raises
    -------
    RuntimeError
        There's no log entry with specified `chopper_speed`
    """
    frequencies = workspace.getRun().getLogData(pulse_freq)
    av_freq = frequencies.getStatistics().mean
    return 1.0e6 / av_freq


def frame_skipping(workspace, chopper_speed='Speed1', pulse_freq='frequency'):
    """
    Verify if neutrons were collected in frame-skipping mode

    Parameters
    ----------
    workspace: Mantid.Workspace
    chopper_speed: str
        Entry name in the log containing the chopper angular speed
    pulse_freq: str
        Entry name in the log containing the pulses generation frequencies

    Returns
    -------
    bool

    Raises
    ------
    RuntimeError
        There's no log entries with specified `chopper_speed` or
        `pulse_freq`
    """
    run_logs = workspace.getRun()
    fr = run_logs.getLogData(pulse_freq).getStatistics().mean
    ch = run_logs.getLogData(chopper_speed).getStatistics().mean
    return abs(ch - fr/2.0) < 1.0
