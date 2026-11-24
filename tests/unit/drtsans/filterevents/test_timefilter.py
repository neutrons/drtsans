import numpy as np
import pytest
from unittest.mock import MagicMock, PropertyMock, patch

from drtsans.filterevents.timefilter import (
    TimeIntervalFilter,
    PeriodicTimeFilter,
    _events_end,
    _slice_lengths,
    _splitter_row,
)


def _make_workspace_group(comments, indexes=None):
    # help function to create a mock workspace group with specified comments for each slice.
    # `indexes` are the FilterEvents suffixes of the slices kept in the group (default 1, 2, ...)
    indexes = indexes if indexes is not None else range(1, len(comments) + 1)
    slices = [
        MagicMock(**{"getComment.return_value": c, "name.return_value": f"output_ws_{k}"})
        for c, k in zip(comments, indexes)
    ]
    group = MagicMock()
    group.getNumberOfEntries.return_value = len(slices)
    group.getItem.side_effect = lambda n: slices[n]
    return group


@pytest.mark.parametrize(
    "name, row",
    [("output_ws_1", 0), ("GPSANS_9166_3", 2), ("BIOSANS_12_slice_10", 9)],
)
def test_splitter_row(name, row):
    assert _splitter_row(name) == row


@pytest.mark.parametrize("name", ["output_ws", "output_ws_", "output_ws_0", "output_ws_a1", "12"])
def test_splitter_row_invalid_name(name):
    with pytest.raises(RuntimeError, match="Cannot find the splitter row"):
        _splitter_row(name)


@pytest.mark.parametrize(
    "time_interval, time_offset, run_duration, expected",
    [
        ([10.0, 20.0], 0.0, 73.0, [10.0, 20.0, 20.0, 23.0]),  # last value repeated, final slice absorbs 3 s
        ([10.0, 20.0], 0.0, 70.0, [10.0, 20.0, 20.0, 20.0]),  # no leftover
        ([10.0, 20.0], 5.0, 73.0, [10.0, 20.0, 38.0]),  # offset: 68 s to partition, final slice 20 <= 38 < 40
        ([10.0, 20.0], 0.0, 45.0, [10.0, 35.0]),  # no copy fits, the list's own last slice absorbs 15 s
        ([10.0, 20.0, 30.0], 0.0, 60.0, [10.0, 20.0, 30.0]),  # exact match
        (10.0, 0.0, 73.0, [10.0] * 6 + [13.0]),  # single number, final slice 10 <= 13 < 20
        (2.0, 1.0, 10.0, [2.0, 2.0, 2.0, 3.0]),  # user documentation, Example 2
        (10.0, 0.0, 70.0, [10.0] * 7),
        ([10.0], 0.0, 73.0, [10.0] * 6 + [13.0]),  # one-item list, same as the number
    ],
)
def test_slice_lengths(time_interval, time_offset, run_duration, expected):
    # Act
    lengths, _ = _slice_lengths(time_interval, time_offset, run_duration)
    # Assert
    assert lengths == pytest.approx(expected)
    assert time_offset + sum(lengths) == pytest.approx(run_duration)


@pytest.mark.parametrize(
    "time_interval, run_duration, level, phrase",
    [
        ([10.0, 20.0], 73.0, "warning", "2 slice(s) of 20 s were added and the final slice lasts 23 s"),
        ([10.0, 20.0], 45.0, "warning", "0 slice(s) of 20 s were added and the final slice lasts 35 s"),
        ([10.0, 20.0, 30.0], 35.0, "error", "covers 60 s from the run start, but the run lasts 35 s"),
        (100.0, 30.0, "error", "covers 100 s from the run start, but the run lasts 30 s"),
        (10.0, 73.0, None, None),  # repeating a single number is not reported
        ([10.0, 20.0, 30.0], 60.0, None, None),  # exact match
    ],
)
def test_slice_lengths_report(time_interval, run_duration, level, phrase):
    # Act
    lengths, report = _slice_lengths(time_interval, 0.0, run_duration)
    # Assert
    if level is None:
        assert report is None
    else:
        assert report[0] == level
        assert phrase in report[1]
    if level == "error":
        # unchanged: the run end cuts the slices short
        assert lengths == (list(time_interval) if isinstance(time_interval, list) else [time_interval])


@pytest.mark.parametrize(
    "time_interval, time_offset, run_duration, expected",
    [
        ([0.1] * 10, 0.0, 1.0, [0.1] * 10),  # sums to 0.9999999999999999 with sum(), exactly 1.0 with fsum
        ([0.1] * 3 + [0.7], 0.3, 2.4, [0.1, 0.1, 0.1, 0.7, 1.1]),  # non-integer values and offset, 0.7 <= 1.1 < 1.4
        ([300.0, 900.0], 0.0, 2100.0000001, [300.0, 900.0, 900.0000001]),  # exact multiple up to rounding
        ([300.0, 900.0], 0.0, 2099.9999999, [300.0, 900.0, 899.9999999]),  # just below, not a 2x final slice
        ([300.0, 900.0], 0.0, 2999.5, [300.0, 900.0, 1799.5]),  # half a second short of a copy: absorbed
    ],
)
def test_slice_lengths_rounding(time_interval, time_offset, run_duration, expected):
    # Act
    lengths, _ = _slice_lengths(time_interval, time_offset, run_duration)
    # Assert
    assert lengths == pytest.approx(expected, abs=1e-9)
    assert all(length > 0 for length in lengths)


def _event_workspace(run_start_ns, last_pulse_s, pulse_period_s=1 / 60, pulses=3):
    # mock event workspace whose last pulse is `last_pulse_s` after the run start, with a proton_charge log
    # holding `pulses` entries spaced by `pulse_period_s`
    workspace = MagicMock()
    run = workspace.run.return_value
    run.startTime.return_value.totalNanoseconds.return_value = run_start_ns
    workspace.getPulseTimeMax.return_value.totalNanoseconds.return_value = run_start_ns + int(last_pulse_s * 1e9)
    times = np.datetime64("2020-01-01T00:00:00", "ns") + np.arange(pulses) * np.timedelta64(int(pulse_period_s * 1e9))
    run.getProperty.return_value.times = times
    return workspace


@pytest.mark.parametrize(
    "last_pulse_s, expected",
    [
        (300.048, 300.048 + 1 / 60),  # events of the last pulse arrive up to one pulse period after it
        (299.0, 300.048),  # last pulse well before the run end
        (-5.4e9, 300.048),  # no events: Mantid returns a placeholder pulse time in 1843
    ],
)
@patch("drtsans.filterevents.timefilter.workspace_handle")
def test_events_end(mock_workspace_handle, last_pulse_s, expected):
    # Arrange
    mock_workspace_handle.return_value = _event_workspace(1_000_000_000, last_pulse_s)
    # Act
    events_end = _events_end("workspace", run_duration=300.048)
    # Assert
    assert events_end == pytest.approx(expected, abs=1e-6)


def _workspace_without_proton_charge():
    workspace = _event_workspace(1_000_000_000, 300.048)
    workspace.run.return_value.getProperty.side_effect = RuntimeError("Unknown property search object proton_charge")
    return workspace


@pytest.mark.parametrize(
    "workspace",
    [
        MagicMock(spec=["run"]),  # not an event workspace: no getPulseTimeMax
        _workspace_without_proton_charge(),
        _event_workspace(1_000_000_000, 300.048, pulses=1),  # a single pulse: no pulse period
    ],
    ids=["not an event workspace", "no proton_charge log", "single pulse"],
)
@patch("drtsans.filterevents.timefilter.workspace_handle")
def test_events_end_fallback(mock_workspace_handle, workspace):
    # Arrange
    mock_workspace_handle.return_value = workspace
    # Act and Assert: the run end is used
    assert _events_end("workspace", run_duration=300.048) == 300.048


@patch("drtsans.filterevents.timefilter._events_end", side_effect=lambda _workspace, run_duration: run_duration)
@patch("drtsans.filterevents.timefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.GenerateEventsFilter")
def test_time_interval_filter_generate_filter(mock_gef, mock_samplelogs_cls, mock_events_end):
    # Arrange
    mock_samplelogs_cls.return_value.run_duration = 73.0

    # Act: default offset
    TimeIntervalFilter(MagicMock(), time_interval=10.0).generate_filter()
    # Assert
    _, kwargs = mock_gef.call_args
    assert kwargs["TimeInterval"] == pytest.approx([10.0] * 6 + [13.0])
    assert kwargs["UnitOfTime"] == "Seconds"
    assert kwargs["StartTime"] == "0.0"
    assert float(kwargs["StopTime"]) == pytest.approx(73.0)

    # Act: explicit offset and a list
    TimeIntervalFilter(MagicMock(), time_interval=[10.0, 20.0], time_offset=5.0).generate_filter()
    # Assert
    _, kwargs = mock_gef.call_args
    assert kwargs["StartTime"] == "5.0"
    assert kwargs["TimeInterval"] == pytest.approx([10.0, 20.0, 38.0])
    assert float(kwargs["StopTime"]) == pytest.approx(73.0)


@patch("drtsans.filterevents.timefilter._events_end", return_value=73.02)
@patch("drtsans.filterevents.timefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.GenerateEventsFilter")
def test_time_interval_filter_generate_filter_last_pulse(mock_gef, mock_samplelogs_cls, mock_events_end):
    # Arrange: the events of the last pulse arrive until 73.02 s, after the run end at 73 s
    mock_samplelogs_cls.return_value.run_duration = 73.0
    # Act
    TimeIntervalFilter(MagicMock(), time_interval=[10.0, 20.0]).generate_filter()
    # Assert: the partition uses the run end, and only the final slice reaches past it, to the last event
    _, kwargs = mock_gef.call_args
    assert kwargs["TimeInterval"] == pytest.approx([10.0, 20.0, 20.0, 23.02])
    assert float(kwargs["StopTime"]) == pytest.approx(73.02)


@patch("drtsans.filterevents.timefilter._events_end", return_value=35.02)
@patch("drtsans.filterevents.timefilter.logger")
@patch("drtsans.filterevents.timefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.GenerateEventsFilter")
def test_time_interval_filter_generate_filter_run_shorter(mock_gef, mock_samplelogs_cls, mock_logger, mock_events_end):
    # Arrange
    mock_samplelogs_cls.return_value.run_duration = 35.0
    # Act
    TimeIntervalFilter(MagicMock(), time_interval=[10.0, 20.0, 30.0]).generate_filter()
    # Assert: the list is passed unchanged and Mantid stops after the last event
    _, kwargs = mock_gef.call_args
    assert kwargs["TimeInterval"] == [10.0, 20.0, 30.0]
    assert float(kwargs["StopTime"]) == pytest.approx(35.02)
    mock_logger.error.assert_called_once()


@patch("drtsans.filterevents.timefilter.logger")
@patch("drtsans.filterevents.timefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.GenerateEventsFilter")
def test_time_interval_filter_generate_filter_no_run_duration(mock_gef, mock_samplelogs_cls, mock_logger):
    # Arrange
    type(mock_samplelogs_cls.return_value).run_duration = PropertyMock(side_effect=RuntimeError("no run_end"))
    # Act
    TimeIntervalFilter(MagicMock(), time_interval=[10.0, 20.0]).generate_filter()
    # Assert: the interval goes to Mantid unchanged, without StopTime
    _, kwargs = mock_gef.call_args
    assert kwargs["TimeInterval"] == [10.0, 20.0]
    assert "StopTime" not in kwargs
    mock_logger.warning.assert_called_once()


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.workspace_handle")
@patch("drtsans.filterevents.basefilter.workspace_handle")
@patch("drtsans.filterevents.timefilter.mtd")
def test_time_interval_inject_metadata(
    mock_timefilter_mtd, mock_workspace_handle, mock_timefilter_workspace_handle, mock_samplelogs_cls
):
    run_start_ns = 1_000_000_000
    comments = ["slice_0", "slice_1"]
    mock_workspace_handle.return_value = _make_workspace_group(comments)
    mock_timefilter_workspace_handle.return_value = mock_workspace_handle.return_value

    samplelogs_instances = []
    for _ in comments:
        sl = MagicMock()
        sl.startTime.return_value.totalNanoseconds.return_value = run_start_ns
        samplelogs_instances.append(sl)
    mock_samplelogs_cls.side_effect = samplelogs_instances

    # Splitter table: row n, col 0 = start ns, col 1 = end ns, slices of different lengths
    # Row 0: [run_start + 1.5s, run_start + 11.5s]
    # Row 1: [run_start + 11.5s, run_start + 34.5s]
    def splitter_cell(row, col):
        starts = [run_start_ns + 1.5e9, run_start_ns + 11.5e9]
        ends = [run_start_ns + 11.5e9, run_start_ns + 34.5e9]
        return starts[row] if col == 0 else ends[row]

    splitter_table = MagicMock()
    splitter_table.cell.side_effect = splitter_cell
    mock_timefilter_mtd.__getitem__.return_value = splitter_table

    TimeIntervalFilter(MagicMock(), time_interval=[10.0, 23.0]).inject_metadata("output_ws")

    sl0 = samplelogs_instances[0]
    sl0.insert.assert_any_call("slice_parameter", "relative time from start")
    sl0.insert.assert_any_call("slice_interval", pytest.approx(10.0), "seconds")
    sl0.insert.assert_any_call("slice_start", pytest.approx(1.5), "seconds")
    sl0.insert.assert_any_call("slice_end", pytest.approx(11.5), "seconds")

    # each slice records its own length
    sl1 = samplelogs_instances[1]
    sl1.insert.assert_any_call("slice_interval", pytest.approx(23.0), "seconds")
    sl1.insert.assert_any_call("slice_start", pytest.approx(11.5), "seconds")
    sl1.insert.assert_any_call("slice_end", pytest.approx(34.5), "seconds")


@pytest.mark.parametrize(
    "kept_indexes, expected_windows",
    [
        # leading empty slices removed, e.g. by FilterByTimeStart=25
        ([3, 4, 5, 6], [(20, 30), (30, 40), (40, 50), (50, 60)]),
        # interior empty slice removed, e.g. by a beam-off period
        ([1, 2, 4, 5], [(0, 10), (10, 20), (30, 40), (40, 50)]),
    ],
    ids=["leading_empty", "interior_empty"],
)
@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.timefilter.workspace_handle")
@patch("drtsans.filterevents.basefilter.workspace_handle")
@patch("drtsans.filterevents.timefilter.mtd")
def test_time_interval_inject_metadata_after_empty_slices_removed(
    mock_timefilter_mtd,
    mock_workspace_handle,
    mock_timefilter_workspace_handle,
    mock_samplelogs_cls,
    kept_indexes,
    expected_windows,
):
    # Arrange: 60 s run sliced every 10 s (six splitter rows); only the slices in `kept_indexes` remain in the group
    run_start_ns = 1_000_000_000
    group = _make_workspace_group([f"slice_{k}" for k in kept_indexes], indexes=kept_indexes)
    mock_workspace_handle.return_value = group
    mock_timefilter_workspace_handle.return_value = group

    samplelogs_instances = []
    for _ in kept_indexes:
        sl = MagicMock()
        sl.startTime.return_value.totalNanoseconds.return_value = run_start_ns
        samplelogs_instances.append(sl)
    mock_samplelogs_cls.side_effect = samplelogs_instances

    splitter_table = MagicMock()
    splitter_table.cell.side_effect = lambda row, col: run_start_ns + (row + col) * 10e9
    mock_timefilter_mtd.__getitem__.return_value = splitter_table

    # Act
    TimeIntervalFilter(MagicMock(), time_interval=10.0).inject_metadata("output_ws")

    # Assert: each slice reports the window of its own splitter row, not the row at its position in the group
    for sl, (start, end) in zip(samplelogs_instances, expected_windows):
        sl.insert.assert_any_call("slice_start", pytest.approx(start), "seconds")
        sl.insert.assert_any_call("slice_end", pytest.approx(end), "seconds")


# ------------------
# PeriodicTimeFilter
# ------------------


@patch("drtsans.filterevents.logfilter.GenerateEventsFilter")
@patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log")
def test_periodic_time_filter_generate_filter(mock_create_log, mock_gef):
    PeriodicTimeFilter(MagicMock(), time_interval=5.0, time_period=60.0).generate_filter()
    mock_gef.assert_called_once()
    _, kwargs = mock_gef.call_args
    assert kwargs["LogName"] == "periodic_time_slicing"
    assert kwargs["LogValueInterval"] == 1


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.workspace_handle")
@patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log")
def test_periodic_time_inject_metadata(mock_create_log, mock_workspace_handle, mock_samplelogs_cls):
    # Comment format parsed by the regex: "...From.X.To.Y.Value..."
    comments = [
        "Splitter.From.0.To.1.Value.0",
        "Splitter.From.1.To.2.Value.1",
    ]
    mock_workspace_handle.return_value = _make_workspace_group(comments)
    samplelogs_instances = [MagicMock() for _ in comments]
    mock_samplelogs_cls.side_effect = samplelogs_instances

    PeriodicTimeFilter(MagicMock(), time_interval=5.0, time_period=60.0).inject_metadata("output_ws")

    sl0 = samplelogs_instances[0]
    sl0.insert.assert_any_call("slice_parameter", "periodic time from start")
    sl0.insert.assert_any_call("slice_interval", 5.0)
    sl0.insert.assert_any_call("slice_period", 60.0)
    # log_start=0, log_end=1  →  0*5=0 s, 1*5=5 s
    sl0.insert.assert_any_call("slice_start", 0.0, "seconds")
    sl0.insert.assert_any_call("slice_end", 5.0, "seconds")

    sl1 = samplelogs_instances[1]
    # log_start=1, log_end=2  →  1*5=5 s, 2*5=10 s
    sl1.insert.assert_any_call("slice_start", 5.0, "seconds")
    sl1.insert.assert_any_call("slice_end", 10.0, "seconds")
