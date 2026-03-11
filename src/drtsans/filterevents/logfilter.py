"""
Sample log-based event filtering strategy.

This module provides filtering based on sample log values, allowing events to be
grouped by changes in experimental parameters recorded in the sample logs.
"""

import re

from drtsans.filterevents.basefilter import FilterStrategy


class LogValueFilter(FilterStrategy):
    """
    Filter events based on sample log values.

    This strategy splits event data based on changes in a sample log parameter.
    Events are grouped into slices where the log value changes by a specified interval.
    This is useful for analyzing data as a function of experimental conditions
    (e.g., temperature, magnetic field, sample position).

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered
    log_name : str
        Name of the sample log to use for filtering (e.g., 'ProtonCharge', 'SampleTemp')
    log_value_interval : float
        Interval of log value changes that define slice boundaries

    Attributes
    ----------
    log_name : str
        The sample log used for filtering
    log_value_interval : float
        The log value interval used for slicing

    Notes
    -----
    The log must exist in the workspace and must be a time series property.
    Events are grouped into intervals [min_value, min_value + interval),
    [min_value + interval, min_value + 2*interval), etc.

    If the log has units, they are preserved in the metadata.
    """

    def __init__(self, workspace, log_name: str, log_value_interval: float):
        """
        Initialize the log value filter.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to filter
        log_name : str
            Name of the sample log to use for filtering
        log_value_interval : float
            Interval of log value changes for slicing
        """
        super().__init__(workspace)
        self.log_name = log_name
        self.log_value_interval = log_value_interval

    def generate_filter(self) -> dict:
        """
        Generate filter parameters for log-based splitting.

        Returns
        -------
        dict
            Parameters for GenerateEventsFilter including:
            - 'LogName': Name of the sample log
            - 'LogValueInterval': The log value interval
        """
        return {"LogName": self.log_name, "LogValueInterval": self.log_value_interval}

    def inject_metadata(self, output_workspace: str) -> None:
        """
        Inject metadata into all log-filtered slices.

        Adds common metadata (slice number, total slices) and log-specific metadata
        (log parameter name, interval, start/end values with units) to each workspace
        in the group.

        Parameters
        ----------
        output_workspace : str
            Name of the workspace group containing the filtered workspaces
        """
        for _, samplelogs, slice_info in self._inject_common_metadata(output_workspace):
            samplelogs.insert("slice_parameter", self.log_name)
            samplelogs.insert("slice_interval", self.log_value_interval)

            slice_start, slice_end = re.sub(r".*\.From\.|\.Value.*", "", slice_info).split(".To.")

            log_units = samplelogs[self.log_name].units if self.log_name in samplelogs else ""

            if log_units:
                samplelogs.insert("slice_start", float(slice_start), log_units)
                samplelogs.insert("slice_end", float(slice_end), log_units)
            else:
                samplelogs.insert("slice_start", float(slice_start))
                samplelogs.insert("slice_end", float(slice_end))
