"""
Time-based event filtering strategies.

This module provides filtering strategies that slice event data based on time intervals,
including both regular and periodic time slicing.
"""

import re
from typing import Union, List
from mantid.simpleapi import mtd

from drtsans.filterevents.basefilter import FilterStrategy
from drtsans.filterevents.logfilter import LogValueFilter
from drtsans.samplelogs import SampleLogs


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

    def generate_filter(self) -> dict:
        """
        Generate filter parameters for time-based splitting.

        Returns
        -------
        dict
            Parameters for GenerateEventsFilter including:
            - 'StartTime': String representation of time offset
            - 'TimeInterval': The time interval(s)
            - 'UnitOfTime': Always 'Seconds' for time filtering
        """
        return {"StartTime": str(self.time_offset), "TimeInterval": self.time_interval, "UnitOfTime": "Seconds"}

    def inject_metadata(self, output_workspace: str) -> None:
        """
        Inject metadata into all time-filtered slices.

        Adds common metadata (slice number, total slices) and time-specific metadata
        (relative times with units) to each workspace in the group.

        Parameters
        ----------
        output_workspace : str
            Name of the workspace group containing the filtered workspaces
        """
        for n, samplelogs, _ in self._inject_common_metadata(output_workspace):
            samplelogs.insert("slice_parameter", "relative time from start")
            samplelogs.insert("slice_interval", self.time_interval)

            start_time_ns = mtd[self.splitter_workspace].cell(n, 0)
            end_time_ns = mtd[self.splitter_workspace].cell(n, 1)
            run_start_ns = samplelogs.startTime().totalNanoseconds()

            samplelogs.insert("slice_start", (start_time_ns - run_start_ns) / 1e9, "seconds")
            samplelogs.insert("slice_end", (end_time_ns - run_start_ns) / 1e9, "seconds")


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

    def generate_filter(self) -> dict:
        """
        Generate filter parameters for periodic time slicing.

        Delegates to the internal LogValueFilter to use log-based filtering
        with the synthetic periodic log.

        Returns
        -------
        dict
            Parameters for GenerateEventsFilter including:
            - 'LogName': Name of the synthetic periodic log
            - 'LogValueInterval': Always 1 (group by log value)
        """
        return self._log_filter.generate_filter()

    def inject_metadata(self, output_workspace: str) -> None:
        """
        Inject metadata into all periodically-filtered slices.

        Adds common metadata and periodic time-specific information to each
        workspace in the group.

        Parameters
        ----------
        output_workspace : str
            Name of the workspace group containing the filtered workspaces
        """
        for _, samplelogs, slice_info in self._inject_common_metadata(output_workspace):
            slice_start_str, slice_end_str = re.sub(r".*\.From\.|\.Value.*", "", slice_info).split(".To.")
            log_start = float(slice_start_str)
            log_end = float(slice_end_str)

            samplelogs.insert("slice_parameter", "periodic time from start")
            samplelogs.insert("slice_interval", self.time_interval)
            samplelogs.insert("slice_period", self.time_period)
            samplelogs.insert("slice_start", log_start * self.time_interval, "seconds")
            samplelogs.insert("slice_end", log_end * self.time_interval, "seconds")
