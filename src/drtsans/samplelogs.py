# package imports

# third party imports
from mantid.api import Run, MatrixWorkspace
from mantid.kernel import (
    BoolTimeSeriesProperty,
    DateAndTime,
    FloatTimeSeriesProperty,
    Int64TimeSeriesProperty,
    StringTimeSeriesProperty,
)
from mantid.simpleapi import mtd
import numpy as np

# standard library imports
from typing import List, Union

TimeSeriesProperty = Union[
    BoolTimeSeriesProperty, FloatTimeSeriesProperty, Int64TimeSeriesProperty, StringTimeSeriesProperty
]
SECONDS_TO_NANOSECONDS = 1.0e09  # from seconds to nanoseconds


def time_series(
    name: str,
    elapsed_times: Union[np.ndarray, List[DateAndTime]],
    values: List[Union[bool, float, int, str]],
    start_time: Union[DateAndTime, str] = "2000-01-01T00:00:00",
    unit: str = "",
) -> TimeSeriesProperty:
    r"""
    Create time series log

    Parameters
    ----------
    name
        log entry name
    elapsed_times
        List of elapsed times after ```start_time```, in seconds.
    start_time
        Starting time for the run
    values
        List of log values, same length as the list of times
    unit
        Log unit

    Returns
    -------
    mantid log that can be added as a property to a Run object

    Raises
    ------
    ValueError
        If items in `sequence `values`` is not one of ``bool, ``float``, ``int``, ``str``
    """
    if len(elapsed_times) != len(values):
        raise ValueError("elapsed_times and values must have the same length")

    python_types = (bool, float, int, str, object)
    mantid_types = (
        BoolTimeSeriesProperty,
        FloatTimeSeriesProperty,
        Int64TimeSeriesProperty,
        StringTimeSeriesProperty,
        FloatTimeSeriesProperty,
    )
    first_value = values[0]
    for python_type, mantid_type in zip(python_types, mantid_types):
        if isinstance(first_value, python_type):
            series_property = mantid_type(name)
            break
    else:
        raise ValueError(f"Cannot create time series log for values of type {type(first_value)}")

    start = start_time
    if isinstance(start_time, str):
        start = DateAndTime(start_time)

    # Insert one pair of (time, elapsed_time) at a time
    for elapsed_time, value in zip(elapsed_times, values):
        series_property.addValue(start + int(elapsed_time * SECONDS_TO_NANOSECONDS), value)

    series_property.units = unit
    return series_property


def periodic_index_log(
    period: float,
    interval: float,
    duration: float,
    run_start: Union[str, DateAndTime],
    offset: float = 0.0,
    step: float = 1.0,
    name: str = "periodic_index",
) -> FloatTimeSeriesProperty:
    r"""
    Generate a periodic log whose values are integers ranging from 0 to ``ceil(period / interval)``
    using timeSliceXXX values from the reduction configuration.

    The first log entry is at ``run_start + offset`` with value 0. The next entry at
    ``run_start + offset + interval`` with value 1, and so on. The log wraps around
    at ``run_start + offset + period`` so the next value is again 0.

    Parameters
    ----------
    period
        Period of the log, in seconds
    interval
        Interval between consecutive log entries, in seconds
    duration
        Duration of the log from ``run_start``, in seconds
    run_start
        Time of the first log entry, unless ``offset`` is also specified
    offset
        Time of the first log entry after ``run_start``, in seconds
    step
        Absolute value of the change in the log value between two consecutive entries
    name
        Name of the log

    Returns
    -------
    A Mantid timeseries property object which can be attached to a Run object.

    """

    # times at which we insert a new log entry
    times = np.arange(offset, duration, interval)

    # number of periods in the duration
    # plus 1 for any remainder (via ceil)
    period_count = int(np.ceil((duration - offset) / period))

    # array of values in each period, scaled by the step
    values_in_period = step * np.arange(0, np.ceil(period / interval))

    # repeat the values in a period up to the number of periods,
    # then truncate to the length of times, then cast to list
    entries = np.tile(values_in_period, period_count)[: len(times)].tolist()

    assert len(times) == len(entries), (
        f"times and entries must have the same length: len(times) {len(times)} != len(entries) {len(entries)}"
    )

    return time_series(name, times, entries, run_start, unit="")


class SampleLogs(object):
    r"""
    Log reader, a bit more pythonic

    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name, run number
    """

    def __init__(self, source):
        self._ws = None
        self._run = self.find_run(source)

    def __getitem__(self, item):
        if item in self._run.keys():
            return self._run[item]
        raise KeyError('"{}" not found in sample logs'.format(item))

    def __getattr__(self, item):
        _run = self.__dict__["_run"]
        try:
            return getattr(_run, item)
        except AttributeError:
            if item in _run.keys():
                return _run.getProperty(item)
            else:
                raise AttributeError('"{}" not found in sample logs'.format(item))

    def __contains__(self, item):
        """Called when using python's ``in`` operation"""
        return item in self._run

    def insert(self, name: str, value: Union[str, int, float, list, TimeSeriesProperty], unit: str = None):
        r"""
        Wrapper to Mantid AddSampleLog algorithm

        The following properties of AddSampleLog are determined by
        inspection of `value`: LogText, LogType, NumberType

        Parameters
        ----------
        name
            log entry name. If ``value`` is one of Mantid's time series property objects, it is
            expected that this is the name of the series.
        value
            Value to insert. If ``value`` is a ``list`` object, it will insert only the first value of the list.
        unit: str
            Log unit. If ``value`` is one of Mantid's time series property objects, it is
            expected that this is the unit of the series.
        """
        if not unit:
            unit = ""
        if isinstance(value, list):
            value = value[0]  # copies AddSampleLog behavior
        if isinstance(
            value, (BoolTimeSeriesProperty, FloatTimeSeriesProperty, Int64TimeSeriesProperty, StringTimeSeriesProperty)
        ):
            assert name == value.name
            assert unit == value.units

        self._ws.mutableRun().addProperty(name, value, unit, replace=True)

    def insert_time_series(self, name, elapsed_times, values, start_time="2000-01-01T00:00:00", unit=""):
        r"""
        Insert a ~mantid.kernel.FloatTimeSeriesProperty in the logs

        Parameters
        ----------
        name: str
            log entry name
        elapsed_times: list
            List of elapsed times after ```start_time```, in seconds.
        values: list
            List of log values, same length as the list of times
        start_time: str
            Starting time for the run
        unit: str
            units of the log values
        """
        series_property = time_series(name, elapsed_times, values, start_time, unit)
        self._ws.mutableRun().addProperty(name, series_property, units=unit, replace=True)

    @property
    def mantid_logs(self):
        return self._run

    @property
    def workspace(self):
        return self._ws

    def single_value(self, log_key, operation=np.mean):
        r"""Cast the log entry to a single value

        Parameters
        ----------
        log_key: str
            the name of the log
        operation: funct
            The casting operation to perform on the values when the log entry is a time series.

        Returns
        -------
        str, float, None
            The single value associated with the log entry.

        Raises
        ------
        RuntimeError
            When the log_key is not found in the sample logs
        """
        _run = self.__dict__["_run"]
        value = _run[log_key].value
        if isinstance(value, str):
            return value
        else:
            return float(operation(value))

    def find_run(self, other):
        r"""
        Retrieve the Run object

        Parameters
        ----------
        other: Run, str, MatrixWorkspace

        Returns
        -------
        Run
            Reference to the run object
        """

        def from_ws(ws):
            self._ws = ws
            return ws.getRun()

        def from_run(a_run):
            return a_run

        def from_string(s):
            # see if it is a file
            if s in mtd:
                return self.find_run(mtd[s])
            else:
                raise RuntimeError("{} is not a valid workspace name".format(s))

        dispatch = {Run: from_run, MatrixWorkspace: from_ws, str: from_string}

        # If others is not None: raise exception
        if other is None:
            a = "[DEBUG 187] dispatch: {}".format(dispatch.items())
            b = "[DEBUG 187] other: {}".format(other)
            raise NotImplementedError("{}\n{}".format(a, b))

        finders = [v for k, v in dispatch.items() if isinstance(other, k)]
        if len(finders) == 0:
            # In case no items found
            raise RuntimeError(
                'Input "other" of value {} is not supported to retrieve Mantid "run" object'.format(other)
            )
        finder = finders[0]
        return finder(other)

    def find_log_with_units(self, log_key, unit=None):
        r"""
        Find a log entry in the logs, and ensure it has the right units

        Parameters
        ----------
        log_key: string
                 key of the log to find
        unit: None or string
               units string to enforce

        Returns
        -------
            log value
        """
        if log_key in self.keys():
            if bool(unit) and not self[log_key].units == unit:
                error_msg = "Found %s with wrong units" % log_key
                error_msg += " [%s]" % self[log_key].units
                raise RuntimeError(error_msg)
            return np.average(self[log_key].value)
        raise RuntimeError(f"Could not find {log_key} with unit {unit} in logs: {self.keys()}")

    def records_json(self, log=None, logs=None, default=None):
        r"""Save the `single_value` of the requested log entries to a dictionary, useful for reporting.

        If both `log` and `logs` are passed, the `log` entry is appended to the `logs` list

        Parameters
        ----------
        log: str
            log entry name
        logs: List[str]
            list of log entry names
        default
            Value assigned to log entry if no entry is found

        Returns
        -------
        dict
            log names and single values for each of the requested log entries
        """
        if logs is None:
            logs = list()
        if log is not None:
            logs.append(log)
        output = {}
        for log_entry in logs:
            try:
                value = self.single_value(log_entry)
            except (KeyError, RuntimeError):
                value = default
            output[log_entry] = value
        return output
