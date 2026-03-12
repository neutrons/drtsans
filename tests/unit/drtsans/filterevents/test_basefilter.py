import pytest
from unittest.mock import MagicMock, patch

from drtsans.filterevents import resolve_slicing
from drtsans.filterevents.basefilter import FilterStrategy, create_filter_strategy
from drtsans.filterevents.logfilter import LogValueFilter
from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter


class _ConcreteFilter(FilterStrategy):
    """Minimal concrete subclass used to exercise FilterStrategy behaviour."""

    def generate_filter(self):
        return {"TimeInterval": 10.0, "UnitOfTime": "Seconds"}

    def inject_metadata(self, workspace) -> None:
        for _, samplelogs, slice_info in self._inject_common_metadata(workspace):
            samplelogs.insert("slice_parameter", "test")


def _make_workspace_group(comments):
    # Create a mock workspace group with the specified comments for each slice
    slices = [MagicMock(**{"getComment.return_value": c}) for c in comments]
    group = MagicMock()
    group.getNumberOfEntries.return_value = len(slices)
    group.getItem.side_effect = lambda n: slices[n]
    return group


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

    # no params → ValueError
    with pytest.raises(ValueError, match="No valid filtering parameters"):
        create_filter_strategy(ws)

    # polarized sample → NotImplementedError
    with pytest.raises(NotImplementedError, match="Spin filtering"):
        create_filter_strategy(
            ws,
            reduction_config={
                "sample": {"runNumber": "123"},
                "configuration": {"polarization": {"level": "half"}},
            },
            time_interval=10.0,
        )


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.workspace_handle")
def test_inject_common_metadata_yields_correct_count(mock_workspace_handle, mock_samplelogs_cls):
    comments = ["info_0", "info_1", "info_2"]
    mock_workspace_handle.return_value = _make_workspace_group(comments)
    mock_samplelogs_cls.side_effect = [MagicMock() for _ in comments]

    results = list(_ConcreteFilter(MagicMock())._inject_common_metadata("output_ws"))
    assert len(results) == 3


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.workspace_handle")
def test_inject_common_metadata_inserts_slice_fields(mock_workspace_handle, mock_samplelogs_cls):
    comments = ["info_0", "info_1", "info_2"]
    mock_workspace_handle.return_value = _make_workspace_group(comments)
    mock_samplelogs_cls.side_effect = [MagicMock() for _ in comments]

    strategy = _ConcreteFilter(MagicMock())
    for n, samplelogs, slice_info in strategy._inject_common_metadata("output_ws"):
        assert slice_info == comments[n]
        samplelogs.insert.assert_any_call("slice", n + 1)
        samplelogs.insert.assert_any_call("number_of_slices", len(comments))
        samplelogs.insert.assert_any_call("slice_info", comments[n])


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.workspace_handle")
def test_inject_common_metadata_n_matches_inserted_slice_number(mock_workspace_handle, mock_samplelogs_cls):
    """Yielded n must always equal the value passed to insert('slice', n+1)."""
    comments = ["a", "b"]
    mock_workspace_handle.return_value = _make_workspace_group(comments)
    mock_samplelogs_cls.side_effect = [MagicMock(), MagicMock()]

    strategy = _ConcreteFilter(MagicMock())
    for n, samplelogs, _ in strategy._inject_common_metadata("output_ws"):
        inserted_slice = next(call.args[1] for call in samplelogs.insert.call_args_list if call.args[0] == "slice")
        assert inserted_slice == n + 1
