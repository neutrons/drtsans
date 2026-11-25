from typing import Optional, Union, List, Tuple

from drtsans.filterevents.basefilter import FilterStrategy
from drtsans.filterevents.logfilter import LogValueFilter
from drtsans.filterevents.spinfilter import SpinFilter
from drtsans.filterevents.timefilter import PeriodicTimeFilter, TimeIntervalFilter
from drtsans.polarization import polarized_sample
from drtsans.type_hints import MantidWorkspace


def create_filter_strategy(
    workspace: MantidWorkspace,
    reduction_config: Optional[dict] = None,
    time_interval: Optional[Union[float, List[float]]] = None,
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
        Time interval(s) for time-based splitting, in seconds. Numeric strings are accepted and cast to float
    time_offset : float, optional
        Offset for the start time, in seconds. Default is 0.0. A numeric string is cast to float
    time_period : float, optional
        Period for periodic time slicing, in seconds. If provided with
        time_interval, creates a PeriodicTimeFilter. A numeric string is cast to float
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
        If no valid filtering parameters are provided, or if time_period is combined with a list of more than
        one time_interval value

    Notes
    -----
    The function determines which filter to use based on which parameters are provided:

    - time_interval + time_period → PeriodicTimeFilter
    - time_interval → TimeIntervalFilter
    - log_name + log_value_interval → LogValueFilter
    - reduction_config --> SpinFilter
    """

    # Determine filter type based on provided parameters
    if time_interval is not None:
        # Values from a configuration that was not validated (validate=False) can still be strings, such as "60"
        if isinstance(time_interval, (list, tuple)):
            time_interval = [float(value) for value in time_interval]
        else:
            time_interval = float(time_interval)
        time_offset = 0.0 if time_offset is None else float(time_offset)
        time_period = None if time_period is None else float(time_period)
        # A one-item list of time intervals is the same as a single interval
        if isinstance(time_interval, (list, tuple)) and len(time_interval) == 1:
            time_interval = time_interval[0]
        # Periodic slicing needs a single interval. The same check runs at schema validation
        # (`_validate_requires_single_interval` in redparams.py); this one covers direct calls to `load_and_split`.
        if isinstance(time_interval, (list, tuple)) and time_period is not None:
            raise ValueError(
                f"time_period ({time_period} s) requires a single time_interval value, but time_interval is a "
                f"list of {len(time_interval)} values {list(time_interval)}. Periodic slicing cannot use a list "
                "of interval lengths."
            )
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
    elif bool(reduction_config) and polarized_sample(reduction_config):
        return SpinFilter(workspace)
    else:
        raise ValueError(
            "No valid filtering parameters provided. Must specify either:\n"
            "  - time_interval (with optional time_period for periodic slicing)\n"
            "  - log_name and log_value_interval\n"
            "  - reduction_config with a polarized sample (spin-based slicing)"
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
