import pytest
from unittest.mock import MagicMock, patch

from drtsans.filterevents.factory import resolve_slicing, create_filter_strategy
from drtsans.filterevents.logfilter import LogValueFilter
from drtsans.filterevents.spinfilter import SpinFilter
from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter


def test_resolve_slicing():
    # no slicing, polarized sample
    assert resolve_slicing(
        {
            "sample": {"runNumber": "12345"},
            "configuration": {"useTimeSlice": False, "useLogSlice": False, "polarization": {"level": "half"}},
        }
    ) == (False, False, True)

    # both time and log slicing requested
    with pytest.raises(ValueError, match="Can't do both time and log slicing"):
        resolve_slicing(
            {
                "sample": {"runNumber": "12345"},
                "configuration": {"useTimeSlice": True, "useLogSlice": True, "polarization": {"level": "half"}},
            }
        )

    # slicing on multiple summed runs
    with pytest.raises(ValueError, match="Can't do slicing on summed data sets"):
        resolve_slicing(
            {
                "sample": {"runNumber": "1,2,3"},
                "configuration": {"useTimeSlice": True, "useLogSlice": False, "polarization": {"level": "half"}},
            }
        )

    # time slicing on polarized data
    with pytest.raises(NotImplementedError, match="Time or log slicing on polarized data sets is not implemented yet"):
        resolve_slicing(
            {
                "sample": {"runNumber": "123"},
                "configuration": {"useTimeSlice": True, "useLogSlice": False, "polarization": {"level": "half"}},
            }
        )

    # log slicing on polarized data
    with pytest.raises(NotImplementedError, match="Time or log slicing on polarized data sets is not implemented yet"):
        resolve_slicing(
            {
                "sample": {"runNumber": "123"},
                "configuration": {"useTimeSlice": False, "useLogSlice": True, "polarization": {"level": "half"}},
            }
        )


def test_create_filter_strategy():

    ws = MagicMock()

    # regular time interval
    strategy = create_filter_strategy(ws, time_interval=10.0)
    assert isinstance(strategy, TimeIntervalFilter)
    assert strategy.time_interval == 10.0
    assert strategy.time_offset == 0.0

    # time interval with explicit offset
    strategy = create_filter_strategy(ws, time_interval=5.0, time_offset=2.0)
    assert isinstance(strategy, TimeIntervalFilter)
    assert strategy.time_offset == 2.0

    # periodic time filter
    with patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log"):
        strategy = create_filter_strategy(ws, time_interval=5.0, time_period=60.0)
    assert isinstance(strategy, PeriodicTimeFilter)
    assert strategy.time_interval == 5.0
    assert strategy.time_period == 60.0

    # log value filter
    strategy = create_filter_strategy(ws, log_name="SampleTemp", log_value_interval=5.0)
    assert isinstance(strategy, LogValueFilter)
    assert strategy.log_name == "SampleTemp"
    assert strategy.log_value_interval == 5.0

    # spin filter
    mock_sl = MagicMock()
    mock_sl.__contains__ = MagicMock(return_value=True)
    mock_sl.get.return_value.value = 1  # both polarizer and analyzer are active
    with patch("drtsans.filterevents.spinfilter.SampleLogs", return_value=mock_sl):
        strategy = create_filter_strategy(
            ws,
            reduction_config={
                "sample": {"runNumber": "123"},
                "configuration": {"polarization": {"level": "half"}},
            },
        )
    assert isinstance(strategy, SpinFilter)
    assert strategy._active_polarizer is True
    assert strategy._active_analyzer is True

    # no params → ValueError
    with pytest.raises(ValueError, match="No valid filtering parameters"):
        create_filter_strategy(ws)


def test_create_filter_strategy_interval_list():
    # Arrange
    ws = MagicMock()

    # Act
    strategy = create_filter_strategy(ws, time_interval=[10.0, 20.0])
    # Assert
    assert isinstance(strategy, TimeIntervalFilter)
    assert strategy.time_interval == [10.0, 20.0]

    # Act: a one-item list is the same as a single interval
    strategy = create_filter_strategy(ws, time_interval=[10.0])
    # Assert
    assert isinstance(strategy, TimeIntervalFilter)
    assert strategy.time_interval == 10.0


def test_create_filter_strategy_periodic_single_interval():
    # Arrange
    ws = MagicMock()

    # Act: a one-item list with a period
    with patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log"):
        strategy = create_filter_strategy(ws, time_interval=[5.0], time_period=60.0)
    # Assert
    assert isinstance(strategy, PeriodicTimeFilter)
    assert strategy.time_interval == 5.0

    # Act and Assert: a list of several intervals with a period
    with pytest.raises(ValueError, match=r"requires a single time_interval value.+list of 2 values \[5.0, 10.0\]"):
        create_filter_strategy(ws, time_interval=[5.0, 10.0], time_period=60.0)


def test_create_filter_strategy_casts_strings():
    # Arrange: values from a configuration that was not validated
    ws = MagicMock()

    # Act
    with patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log"):
        periodic = create_filter_strategy(ws, time_interval="60", time_offset="5", time_period="120")
    regular = create_filter_strategy(ws, time_interval=["10", 20], time_offset=None)

    # Assert
    assert isinstance(periodic, PeriodicTimeFilter)
    assert (periodic.time_interval, periodic.time_offset, periodic.time_period) == (60.0, 5.0, 120.0)
    assert isinstance(regular, TimeIntervalFilter)
    assert (regular.time_interval, regular.time_offset) == ([10.0, 20.0], 0.0)
