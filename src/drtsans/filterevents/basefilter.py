"""
Base strategy class and factory for event filtering.

This module defines the abstract base class for all event filtering strategies,
along with utility functions for determining which strategy to use based on
reduction configuration parameters.
"""

from abc import ABC, abstractmethod
from typing import Optional, Union, Tuple
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
    ...     def get_metadata(self, slice_idx):
    ...         return {'slice_parameter': 'time'}
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
    def get_metadata(self, slice_idx: int) -> dict:
        """
        Get metadata specific to this filter type for a given slice.

        This method extracts relevant information about the filtering parameters
        and the specific slice to be added to the workspace's sample logs.

        Parameters
        ----------
        slice_idx : int
            Index of the filtered workspace slice (0-based)

        Returns
        -------
        dict
            Dictionary of metadata key-value pairs to add to sample logs.
            Common keys include:
            - 'slice_parameter': Name or description of the slicing parameter
            - 'slice_interval': The interval value used for slicing
            - 'slice_start': Start value/time of this slice
            - 'slice_end': End value/time of this slice

        Examples
        --------
        >>> def get_metadata(self, slice_idx):
        ...     return {
        ...         'slice_parameter': 'relative time from start',
        ...         'slice_interval': self.time_interval,
        ...         'slice_start': slice_idx * self.time_interval,
        ...         'slice_end': (slice_idx + 1) * self.time_interval
        ...     }
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
    workspace: MantidWorkspace, reduction_config: Optional[dict] = None, **kwargs
) -> FilterStrategy:
    """
    Factory function to create the appropriate filter strategy.

    This function analyzes the input parameters and instantiates the correct
    FilterStrategy subclass. It supports both reduction configuration dictionaries
    and direct parameter specifications.

    Parameters
    ----------
    workspace : MantidWorkspace
        Input workspace to filter (workspace object or name)
    reduction_config : dict, optional
        Reduction configuration dictionary. If provided, the filtering type and
        parameters are extracted from this dictionary using :func:`resolve_slicing`.
        Should contain keys like "configuration" with "useTimeSlice", "useLogSlice",
        "timeSliceInterval", "logSliceName", etc.
    **kwargs
        Direct filter parameters. Supported parameters depend on the filter type:

        - For time filtering: ``time_interval``, ``time_offset``
        - For periodic time filtering: ``time_interval``, ``time_period``, ``time_offset``
        - For log filtering: ``log_name``, ``log_value_interval``
        - For spin filtering: ``pv_polarizer_state``, ``pv_analyzer_state``, etc.

    Returns
    -------
    FilterStrategy
        An instantiated filter strategy ready to use

    Raises
    ------
    ValueError
        If no valid filtering parameters are provided

    Examples
    --------
    Create a filter from reduction config:

    >>> config = {
    ...     "configuration": {
    ...         "useTimeSlice": True,
    ...         "timeSliceInterval": 10.0,
    ...         "timeSliceOffset": 5.0
    ...     },
    ...     "sample": {"runNumber": "12345"}
    ... }
    >>> strategy = create_filter_strategy(workspace, reduction_config=config)
    >>> strategy.apply_filter("output_ws")

    Create a filter with direct parameters:

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

    Notes
    -----
    When both ``reduction_config`` and ``kwargs`` are provided, ``reduction_config``
    takes precedence for determining the filter type, but ``kwargs`` can provide
    additional parameters not present in the config.
    """
    # Import here to avoid circular dependencies
    from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter
    from drtsans.filterevents.logfilter import LogValueFilter
    from drtsans.filterevents.spinfilter import SpinFilter

    # Resolve from reduction config if provided
    if reduction_config:
        time_slice, log_slice, spin_slice = resolve_slicing(reduction_config)

        if spin_slice:
            return SpinFilter(workspace, **kwargs)
        elif time_slice:
            # Extract time parameters from config
            config = reduction_config["configuration"]
            time_interval = config.get("timeSliceInterval")
            time_offset = config.get("timeSliceOffset", 0.0)
            time_period = config.get("timeSlicePeriod")

            # Override with kwargs if provided
            time_interval = kwargs.get("time_interval", time_interval)
            time_offset = kwargs.get("time_offset", time_offset)
            time_period = kwargs.get("time_period", time_period)

            if time_period:
                return PeriodicTimeFilter(workspace, time_interval, time_period, time_offset)
            else:
                return TimeIntervalFilter(workspace, time_interval, time_offset)
        elif log_slice:
            # Extract log parameters from config
            config = reduction_config["configuration"]
            log_name = config.get("logSliceName")
            log_value_interval = config.get("logSliceInterval")

            # Override with kwargs if provided
            log_name = kwargs.get("log_name", log_name)
            log_value_interval = kwargs.get("log_value_interval", log_value_interval)

            return LogValueFilter(workspace, log_name, log_value_interval)

    # Fallback: infer from explicit kwargs
    if "time_period" in kwargs and "time_interval" in kwargs:
        return PeriodicTimeFilter(
            workspace, kwargs["time_interval"], kwargs["time_period"], kwargs.get("time_offset", 0.0)
        )
    elif "time_interval" in kwargs:
        return TimeIntervalFilter(workspace, kwargs["time_interval"], kwargs.get("time_offset", 0.0))
    elif "log_name" in kwargs and "log_value_interval" in kwargs:
        return LogValueFilter(workspace, kwargs["log_name"], kwargs["log_value_interval"])
    elif any(key.startswith("pv_") for key in kwargs) or "check_devices" in kwargs:
        return SpinFilter(workspace, **kwargs)

    raise ValueError("No valid filtering parameters provided")
