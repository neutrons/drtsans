"""
Base strategy class and factory for event filtering.

This module defines the abstract base class for all event filtering strategies,
along with utility functions for determining which strategy to use based on
reduction configuration parameters.
"""

from abc import ABC, abstractmethod
from typing import Optional, Union, Tuple, List
from mantid.api import IEventWorkspace
from mantid.simpleapi import GenerateEventsFilter, FilterEvents

from drtsans.polarization import polarized_sample
from drtsans.type_hints import MantidWorkspace


class FilterStrategy(ABC):
    """
    Abstract base class for event filtering strategies.

    This class defines the interface that all concrete filtering strategies must implement.
    It provides a common framework for generating filters, applying them to workspaces,
    and extracting metadata from the filtered results.

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered

    Attributes
    ----------
    workspace : str or IEventWorkspace
        The input workspace being filtered
    splitter_workspace : str
        Name of the temporary workspace holding the splitter table
    info_workspace : str
        Name of the temporary workspace holding filter information

    Examples
    --------
    Concrete strategies should inherit from this class:

    >>> class MyFilter(FilterStrategy):
    ...     def generate_filter(self):
    ...         return {'TimeInterval': 10.0}
    ...     def inject_metadata(self, output_workspace):
    ...         # Add metadata to each slice
    ...         pass
    """

    def __init__(self, workspace: Union[str, IEventWorkspace]):
        """
        Initialize the filter strategy.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to be filtered
        """
        self.workspace = workspace
        self.splitter_workspace = "_filter"
        self.info_workspace = "_info"

    @abstractmethod
    def generate_filter(self) -> Optional[dict]:
        """
        Generate filter parameters for Mantid's GenerateEventsFilter algorithm.

        This method must be implemented by concrete strategy classes to provide
        the specific parameters needed for their filtering approach.

        Returns
        -------
        dict or None
            Parameters to pass to GenerateEventsFilter algorithm.
            Return None if filtering is not applicable (e.g., no polarization devices found).

        Examples
        --------
        >>> def generate_filter(self):
        ...     return {
        ...         'StartTime': '0.0',
        ...         'TimeInterval': self.time_interval,
        ...         'UnitOfTime': 'Seconds'
        ...     }
        """
        pass

    @abstractmethod
    def inject_metadata(self, output_workspace: str) -> None:
        """
        Inject metadata into all sliced workspaces in the output workspace group.

        This method must be implemented by concrete strategy classes to add
        filter-specific metadata to each workspace slice. It should add both
        common metadata (slice number, total slices) and filter-specific metadata
        (parameter name, interval, start/end values with appropriate units).

        Parameters
        ----------
        output_workspace : str
            Name of the workspace group containing the filtered workspaces

        Notes
        -----
        Each implementation should add the following metadata to each slice:
        - 'slice': Slice number (1-based)
        - 'number_of_slices': Total number of slices
        - 'slice_info': Comment from the workspace
        - 'slice_parameter': Name/description of the slicing parameter
        - 'slice_interval': The interval value used
        - 'slice_start': Start value for this slice (with units if applicable)
        - 'slice_end': End value for this slice (with units if applicable)

        Examples
        --------
        >>> def inject_metadata(self, output_workspace):
        ...     workspace_group = mtd[output_workspace]
        ...     for n in range(workspace_group.getNumberOfEntries()):
        ...         slice_ws = workspace_group.getItem(n)
        ...         samplelogs = SampleLogs(slice_ws)
        ...         # Add common metadata
        ...         samplelogs.insert("slice", n + 1)
        ...         samplelogs.insert("number_of_slices", workspace_group.getNumberOfEntries())
        ...         # Add filter-specific metadata
        ...         samplelogs.insert("slice_parameter", "relative time from start")
        ...         # ... add interval, start/end times with units
        """
        pass

    def apply_filter(self, output_workspace: str) -> None:
        """
        Apply the filtering logic to create split workspaces.

        This method can be overridden by subclasses that need custom filtering
        behavior (e.g., polarization filtering with custom splitter tables).
        The default implementation uses Mantid's standard GenerateEventsFilter
        and FilterEvents algorithms.

        Parameters
        ----------
        output_workspace : str
            Name for the output workspace group that will contain filtered workspaces

        Notes
        -----
        After filtering, empty workspaces should be removed from the group.
        Subclasses overriding this method should maintain similar behavior.

        Examples
        --------
        >>> strategy = TimeIntervalFilter(workspace, time_interval=10.0)
        >>> strategy.apply_filter("filtered_data")
        >>> # Creates workspace group "filtered_data" with time-sliced data
        """
        filter_params = self.generate_filter()

        # Allow subclasses to handle no-filter case
        if filter_params is None:
            return

        GenerateEventsFilter(
            InputWorkspace=self.workspace,
            OutputWorkspace=self.splitter_workspace,
            InformationWorkspace=self.info_workspace,
            **filter_params,
        )

        FilterEvents(
            InputWorkspace=self.workspace,
            OutputWorkspaceBaseName=output_workspace,
            SplitterWorkspace=self.splitter_workspace,
            InformationWorkspace=self.info_workspace,
            FilterByPulseTime=True,
            GroupWorkspaces=True,
            OutputWorkspaceIndexedFrom1=True,
        )


def resolve_slicing(reduction_input: dict) -> Tuple[bool, bool, bool]:
    """
    Determine which filtering strategy to use based on reduction configuration.

    This function analyzes the reduction configuration to determine whether time-based,
    log-based, or spin-based (polarization) filtering should be applied. It also
    validates that the configuration is internally consistent.

    Parameters
    ----------
    reduction_input : dict
        Dictionary of reduction configuration parameters. Must contain a "configuration"
        key with sub-keys "useTimeSlice" and "useLogSlice". Must also contain a "sample"
        key with a "runNumber" sub-key.

    Returns
    -------
    tuple of bool
        A 3-tuple (timeslice, logslice, spinslice) indicating which filtering type
        should be applied:

        - timeslice: True if time-based slicing is requested
        - logslice: True if log-based slicing is requested
        - spinslice: True if polarization-based slicing is requested

    Raises
    ------
    ValueError
        - If both time and log slicing are True (mutually exclusive)
        - If slicing is requested on multiple summed runs
    NotImplementedError
        - If time or log slicing is requested on polarized data

    Examples
    --------
    >>> config = {
    ...     "configuration": {"useTimeSlice": True, "useLogSlice": False},
    ...     "sample": {"runNumber": "12345"}
    ... }
    >>> timeslice, logslice, spinslice = resolve_slicing(config)
    >>> assert timeslice and not logslice and not spinslice

    Notes
    -----
    This function imports `polarized_sample` from `drtsans.polarization` to detect
    if the experiment uses polarization. Polarization filtering takes precedence
    over time and log slicing.
    """

    reduction_config = reduction_input["configuration"]
    timeslice = reduction_config.get("useTimeSlice", False)
    logslice = reduction_config.get("useLogSlice", False)

    if timeslice and logslice:
        raise ValueError("Can't do both time and log slicing")

    sample = reduction_input["sample"]["runNumber"]
    multiple_samples = len(sample.split(",")) > 1
    if (timeslice or logslice) and multiple_samples:
        raise ValueError("Can't do slicing on summed data sets")

    polarized = polarized_sample(reduction_input)
    if (timeslice and polarized) or (logslice and polarized):
        raise NotImplementedError("Time or log slicing on polarized data sets is not implemented yet")

    return timeslice, logslice, polarized


def create_filter_strategy(
    workspace: MantidWorkspace,
    reduction_config: Optional[dict] = None,
    time_interval: Optional[Union[float, List]] = None,
    time_offset: float = 0.0,
    time_period: Optional[float] = None,
    log_name: Optional[str] = None,
    log_value_interval: Optional[float] = None,
) -> FilterStrategy:
    """
    Factory function to create the appropriate filter strategy.

    This function analyzes the input parameters and instantiates the correct
    FilterStrategy subclass based on the provided filtering parameters.

    Parameters
    ----------
    workspace : MantidWorkspace
        Input workspace to filter (workspace object or name)
    reduction_config : dict, optional
        Reduction configuration dictionary. Only used to detect if the sample
        is polarized. All other filter parameters should be passed explicitly.
    time_interval : float or list of float, optional
        Time interval(s) for time-based splitting, in seconds
    time_offset : float, optional
        Offset for the start time, in seconds. Default is 0.0
    time_period : float, optional
        Period for periodic time slicing, in seconds. If provided with
        time_interval, creates a PeriodicTimeFilter
    log_name : str, optional
        Name of the sample log for log-based filtering
    log_value_interval : float, optional
        Interval of log value changes for log-based filtering

    Returns
    -------
    FilterStrategy
        An instantiated filter strategy ready to use

    Raises
    ------
    ValueError
        If no valid filtering parameters are provided
    NotImplementedError
        If polarized sample is detected (spin filtering not yet supported)

    Examples
    --------
    Create a time filter:

    >>> strategy = create_filter_strategy(
    ...     workspace,
    ...     time_interval=10.0,
    ...     time_offset=5.0
    ... )
    >>> strategy.apply_filter("output_ws")

    Create a periodic time filter:

    >>> strategy = create_filter_strategy(
    ...     workspace,
    ...     time_interval=5.0,
    ...     time_period=60.0
    ... )
    >>> strategy.apply_filter("output_ws")

    Create a log filter:

    >>> strategy = create_filter_strategy(
    ...     workspace,
    ...     log_name='SampleTemp',
    ...     log_value_interval=5.0
    ... )
    >>> strategy.apply_filter("output_ws")

    Notes
    -----
    The function determines which filter to use based on which parameters
    are provided:

    - time_interval + time_period → PeriodicTimeFilter
    - time_interval → TimeIntervalFilter
    - log_name + log_value_interval → LogValueFilter

    Polarized samples are not yet supported by this factory function.
    """
    # Import here to avoid circular dependencies
    from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter
    from drtsans.filterevents.logfilter import LogValueFilter

    # Check for polarized sample - not yet supported
    if reduction_config and polarized_sample(reduction_config):
        raise NotImplementedError(
            "Spin filtering for polarized samples is not yet implemented in create_filter_strategy"
        )

    # Determine filter type based on provided parameters
    if time_interval is not None:
        # Check if this is periodic time slicing
        if time_period is not None:
            # Periodic time slicing
            return PeriodicTimeFilter(
                workspace,
                time_interval=time_interval,
                time_period=time_period,
                time_offset=time_offset,
            )
        else:
            # Regular time slicing
            return TimeIntervalFilter(
                workspace,
                time_interval=time_interval,
                time_offset=time_offset,
            )
    elif log_name is not None and log_value_interval is not None:
        # Log-based slicing
        return LogValueFilter(
            workspace,
            log_name=log_name,
            log_value_interval=log_value_interval,
        )

    raise ValueError(
        "No valid filtering parameters provided. Must specify either:\n"
        "  - time_interval (with optional time_period for periodic slicing)\n"
        "  - log_name and log_value_interval"
    )
