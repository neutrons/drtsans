"""
Time-based event filtering strategies.

This module provides filtering strategies that slice event data based on time intervals,
including both regular and periodic time slicing.
"""

import math
import re
from typing import List, Optional, Tuple, Union

import numpy as np
from mantid.kernel import Logger
from mantid.simpleapi import GenerateEventsFilter, mtd

from drtsans.dataobjects import workspace_handle
from drtsans.filterevents.basefilter import FilterStrategy
from drtsans.filterevents.logfilter import LogValueFilter
from drtsans.samplelogs import SampleLogs
from drtsans.type_hints import MantidWorkspace

logger = Logger("drtsans.filterevents.timefilter")

# Absorbs floating-point rounding only, in seconds. It is not a tolerance on the slice lengths.
_ROUNDING_TOLERANCE = 1.0e-6


def _slice_lengths(
    time_interval: Union[float, List[float]], time_offset: float, run_duration: float
) -> Tuple[List[float], Optional[Tuple[str, str]]]:
    """
    Partition the run into consecutive slice lengths, applying the last-slice rule.

    The slices start at ``time_offset`` and use the values of ``time_interval`` in order. When the run lasts longer,
    the last value ``t`` is repeated, and the final slice absorbs the leftover of the run, so that its length ``x``
    satisfies ``t <= x < 2t``. A single number is treated as a one-item list. When the run is shorter than the values,
    the lengths are returned unchanged: the run end cuts the slices short.

    Parameters
    ----------
    time_interval : float or list of float
        Slice length(s), in seconds
    time_offset : float
        Start of the first slice, in seconds from the run start
    run_duration : float
        Duration of the whole run, in seconds, independent of any load-time cutoffs

    Returns
    -------
    tuple
        The list of slice lengths, in seconds, and either ``None`` or a ``(level, message)`` pair to log, with level
        ``"warning"`` when a list was extended and ``"error"`` when the run is shorter than the values.
    """
    is_list = isinstance(time_interval, (list, tuple))
    lengths = [float(value) for value in time_interval] if is_list else [float(time_interval)]
    available = run_duration - time_offset  # time from the start of the first slice to the run end
    covered = math.fsum(lengths)
    if covered > available + _ROUNDING_TOLERANCE:
        message = (
            f"timeSliceInterval {time_interval} covers {time_offset + covered:.6g} s from the run start, but the run "
            f"lasts {run_duration:.6g} s; the run end cuts the last slice short and later slices are not produced"
        )
        return lengths, ("error", message)

    remaining = available - covered
    if remaining <= _ROUNDING_TOLERANCE:
        return lengths, None  # the values cover the run exactly
    last = lengths[-1]
    # copies of the last value that fit in the remaining time; the tolerance keeps a remainder that is a multiple of
    # `last` up to rounding from turning into a final slice of twice its length
    copies = math.floor((remaining + _ROUNDING_TOLERANCE) / last)
    lengths.extend([last] * copies)
    lengths[-1] += remaining - copies * last  # the final slice absorbs the leftover, ending at the run end
    if not (is_list and len(time_interval) > 1):
        return lengths, None  # repeating a single number is its meaning, nothing to report
    message = (
        f"timeSliceInterval {time_interval} covers {time_offset + covered:.6g} s from the run start, but the run "
        f"lasts {run_duration:.6g} s; {copies} slice(s) of {last:.6g} s were added and the final slice lasts "
        f"{lengths[-1]:.6g} s"
    )
    return lengths, ("warning", message)


def _events_end(workspace: MantidWorkspace, run_duration: float) -> float:
    """
    End of the last pulse, in seconds from the run start, and never earlier than the run end.

    ``FilterEvents`` places each event by its pulse time plus its time of flight, so the events of the last pulse
    arrive after that pulse by up to one pulse period. The final slice must reach past them. The largest time of
    flight is not used: raw data can hold a few corrupt events with a time of flight of hundreds of seconds.

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The workspace to be split
    run_duration : float
        Duration of the run, in seconds

    Returns
    -------
    float
        The later of ``run_duration`` and the last pulse time plus one pulse period. The pulse period is the median
        spacing of the ``proton_charge`` log entries. ``run_duration`` when the workspace has no pulse times or no
        such log.
    """
    try:
        events = workspace_handle(workspace)
        run = events.run()
        last_pulse_ns = events.getPulseTimeMax().totalNanoseconds() - run.startTime().totalNanoseconds()
        pulse_times = np.asarray(run.getProperty("proton_charge").times, dtype="datetime64[ns]")
    except (AttributeError, KeyError, RuntimeError):  # not an event workspace, or no proton_charge log
        return run_duration
    if pulse_times.size < 2:
        return run_duration
    pulse_period_ns = float(np.median(np.diff(pulse_times).astype("int64")))
    # an event workspace without events returns a placeholder pulse time long before the run start
    return max(run_duration, (last_pulse_ns + pulse_period_ns) / 1e9)


def _splitter_row(workspace_name: str) -> int:
    """
    Return the splitter-table row whose events a sliced workspace holds.

    ``FilterEvents`` (with ``OutputWorkspaceIndexedFrom1=True``) names the workspace holding the events
    of splitter row ``k - 1`` as ``<base>_k``.

    Parameters
    ----------
    workspace_name : str
        Name of a workspace produced by ``FilterEvents``, e.g. ``GPSANS_9166_3``

    Returns
    -------
    int
        Zero-based row index in the splitter table

    Raises
    ------
    RuntimeError
        If the name does not end with ``_k``, ``k`` being a positive integer
    """
    _, separator, suffix = workspace_name.rpartition("_")
    if not separator or not suffix.isdigit() or int(suffix) < 1:
        raise RuntimeError(
            f"Cannot find the splitter row of sliced workspace '{workspace_name}': "
            "expected a name ending in '_k', k being a positive integer"
        )
    return int(suffix) - 1


class TimeIntervalFilter(FilterStrategy):
    """
    Filter events by regular time intervals.

    This strategy splits event data into time slices of specified duration(s).
    Time intervals can be either a single value (creating equal-duration slices)
    or a list of values (creating variable-duration slices).

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered
    time_interval : float or list of float
        Time interval(s) for splitting, in seconds. If a single float, all slices
        will have the same duration. If a list, each slice can have a different duration.
    time_offset : float, optional
        Offset to add to the start time of the first slice, in seconds. Default is 0.0.

    Attributes
    ----------
    time_interval : float or list of float
        The time interval(s) used for slicing
    time_offset : float
        The time offset applied to the first slice

    Notes
    -----
    The time interval is measured from the start of the run (after applying offset).
    Events outside the defined time windows are discarded.
    """

    def __init__(self, workspace, time_interval: Union[float, List[float]], time_offset: float = 0.0):
        """
        Initialize the time interval filter.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to filter
        time_interval : float or list of float
            Time interval(s) for splitting, in seconds
        time_offset : float, optional
            Offset for the start time of the first slice, in seconds
        """
        super().__init__(workspace)
        self.time_interval = time_interval
        self.time_offset = time_offset

    def generate_filter(self) -> None:
        """
        Generate and apply the time-based event filter.

        Calls Mantid's GenerateEventsFilter directly, writing the splitter
        and information workspaces used by ``apply_filter``. The slice lengths are
        computed here (see ``_slice_lengths``) and passed to Mantid as an explicit list
        ending at ``StopTime``, because Mantid would otherwise start a list over and
        leave a shorter last slice.
        """
        options = {"TimeInterval": self.time_interval}
        try:
            # whole-run duration, from the run start and end logs, which load-time cutoffs do not change
            run_duration = SampleLogs(self.workspace).run_duration
        except RuntimeError:
            logger.warning("Cannot read the run duration; the time slices are left to Mantid's GenerateEventsFilter")
        else:
            lengths, report = _slice_lengths(self.time_interval, self.time_offset, run_duration)
            if report is not None:
                level, message = report
                getattr(logger, level)(message)
            slices_end = self.time_offset + math.fsum(lengths)
            events_end = _events_end(self.workspace, run_duration)
            if slices_end <= run_duration + _ROUNDING_TOLERANCE:
                # the slices end at the run end: the final slice also takes the events of the last pulse
                lengths[-1] += events_end - slices_end
                slices_end = events_end
            stop_time = min(slices_end, events_end)
            options = {"TimeInterval": lengths, "StopTime": str(stop_time)}
        GenerateEventsFilter(
            InputWorkspace=self.workspace,
            OutputWorkspace=self.splitter_workspace,
            InformationWorkspace=self.info_workspace,
            StartTime=str(self.time_offset),
            UnitOfTime="Seconds",
            **options,
        )

    def inject_metadata(self, workspace: MantidWorkspace) -> None:
        """
        Inject metadata into all time-filtered slices.

        Adds common metadata (slice number, total slices) and time-specific metadata
        (relative times with units) to each workspace in the group.

        Parameters
        ----------
        workspace : MantidWorkspace
            The workspace group (or its name) containing the filtered workspaces
        """
        workspace_group = workspace_handle(workspace)
        for n, samplelogs, _ in self._inject_common_metadata(workspace):
            samplelogs.insert("slice_parameter", "relative time from start")

            # Empty slices have already been removed from the group, so the slice's position in the group
            # no longer matches its row in the splitter table. The row is encoded in the workspace name.
            row = _splitter_row(workspace_group.getItem(n).name())
            start_time_ns = mtd[self.splitter_workspace].cell(row, 0)
            end_time_ns = mtd[self.splitter_workspace].cell(row, 1)
            run_start_ns = samplelogs.startTime().totalNanoseconds()

            samplelogs.insert("slice_start", (start_time_ns - run_start_ns) / 1e9, "seconds")
            samplelogs.insert("slice_end", (end_time_ns - run_start_ns) / 1e9, "seconds")
            # each slice records its own length, from the integer nanosecond times of its splitter row
            samplelogs.insert("slice_interval", (end_time_ns - start_time_ns) / 1e9, "seconds")


class PeriodicTimeFilter(FilterStrategy):
    """
    Filter events by periodic time intervals.

    This strategy handles periodic time slicing by creating a synthetic sample log
    that tracks the periodic structure, then delegating to LogValueFilter for the
    actual filtering. Events that occur at the same relative position within different
    periods are grouped together. This is useful for analyzing periodic phenomena or
    for averaging over multiple cycles.

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered
    time_interval : float
        Duration of each slice within a period, in seconds
    time_period : float
        Duration of the full period, in seconds. Must be a multiple of time_interval.
    time_offset : float, optional
        Offset to add to the start time, in seconds. Default is 0.0.

    Attributes
    ----------
    time_interval : float
        Duration of each slice within a period
    time_period : float
        Duration of the full period
    time_offset : float
        Time offset applied
    log_name : str
        Name of the synthetic log created ('periodic_time_slicing')

    Notes
    -----
    This filter creates a temporary sample log that emulates the periodicity.
    Within each period, the log value increments by 1 for each time_interval.
    The log value resets to 0 at the start of each new period.

    Mantid's GenerateEventsFilter doesn't natively support periodic time slicing,
    so this approach converts the problem into log-based slicing by delegating
    to LogValueFilter.
    """

    def __init__(self, workspace, time_interval: float, time_period: float, time_offset: float = 0.0):
        """
        Initialize the periodic time filter.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to filter
        time_interval : float
            Duration of each slice within a period, in seconds
        time_period : float
            Duration of the full period, in seconds
        time_offset : float, optional
            Offset for the start time, in seconds
        """
        super().__init__(workspace)
        self.time_interval = time_interval
        self.time_period = time_period
        self.time_offset = time_offset
        self.log_name = "periodic_time_slicing"

        # Create the synthetic periodic log
        self._create_periodic_log()

        # Create a LogValueFilter to handle the actual filtering
        self._log_filter: LogValueFilter = LogValueFilter(
            workspace=workspace,
            log_name=self.log_name,
            log_value_interval=1,  # Group by each log value change
        )

    def _create_periodic_log(self):
        """
        Create a synthetic sample log that emulates periodic time structure.

        This method inserts a log into the workspace where the value increments
        by 1 every time_interval seconds, resetting to 0 at the start of each
        new period.
        """
        from drtsans.samplelogs import periodic_index_log

        sample_logs = SampleLogs(self.workspace)

        try:
            run_start = sample_logs.run_start.value
        except AttributeError:
            run_start = sample_logs.start_time.value

        log = periodic_index_log(
            period=self.time_period,
            interval=self.time_interval,
            duration=sample_logs.run_duration,
            run_start=run_start,
            offset=self.time_offset,
            step=1.0,
            name=self.log_name,
        )
        sample_logs.insert(name=self.log_name, value=log)

    def generate_filter(self) -> None:
        """
        Generate and apply the periodic time-based event filter.

        Delegates to the internal ``LogValueFilter`` which calls Mantid's
        GenerateEventsFilter directly using the synthetic periodic log.
        """
        self._log_filter.generate_filter()

    def inject_metadata(self, workspace: MantidWorkspace) -> None:
        """
        Inject metadata into all periodically-filtered slices.

        Adds common metadata and periodic time-specific information to each
        workspace in the group.

        Parameters
        ----------
        workspace : MantidWorkspace
            The workspace group (or its name) containing the filtered workspaces
        """
        for _, samplelogs, slice_info in self._inject_common_metadata(workspace):
            slice_start_str, slice_end_str = re.sub(r".*\.From\.|\.Value.*", "", slice_info).split(".To.")
            log_start = float(slice_start_str)
            log_end = float(slice_end_str)

            samplelogs.insert("slice_parameter", "periodic time from start")
            samplelogs.insert("slice_interval", self.time_interval)
            samplelogs.insert("slice_period", self.time_period)
            samplelogs.insert("slice_start", log_start * self.time_interval, "seconds")
            samplelogs.insert("slice_end", log_end * self.time_interval, "seconds")
