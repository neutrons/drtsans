import pytest
from unittest.mock import MagicMock, patch

from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter, _splitter_row


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


@patch("drtsans.filterevents.timefilter.GenerateEventsFilter")
def test_time_interval_filter_generate_filter(mock_gef):
    # default offset
    TimeIntervalFilter(MagicMock(), time_interval=10.0).generate_filter()
    _, kwargs = mock_gef.call_args
    assert kwargs["TimeInterval"] == 10.0
    assert kwargs["UnitOfTime"] == "Seconds"
    assert kwargs["StartTime"] == "0.0"

    # explicit offset
    TimeIntervalFilter(MagicMock(), time_interval=5.0, time_offset=3.0).generate_filter()
    _, kwargs = mock_gef.call_args
    assert kwargs["StartTime"] == "3.0"
    assert kwargs["TimeInterval"] == 5.0


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

    # Splitter table: row n, col 0 = start ns, col 1 = end ns
    # Row 0: [run_start + 1.5s, run_start + 11.5s]
    # Row 1: [run_start + 11.5s, run_start + 21.5s]
    def splitter_cell(row, col):
        starts = [run_start_ns + 1.5e9, run_start_ns + 11.5e9]
        ends = [run_start_ns + 11.5e9, run_start_ns + 21.5e9]
        return starts[row] if col == 0 else ends[row]

    splitter_table = MagicMock()
    splitter_table.cell.side_effect = splitter_cell
    mock_timefilter_mtd.__getitem__.return_value = splitter_table

    TimeIntervalFilter(MagicMock(), time_interval=10.0).inject_metadata("output_ws")

    sl0 = samplelogs_instances[0]
    sl0.insert.assert_any_call("slice_parameter", "relative time from start")
    sl0.insert.assert_any_call("slice_interval", 10.0)
    sl0.insert.assert_any_call("slice_start", pytest.approx(1.5), "seconds")
    sl0.insert.assert_any_call("slice_end", pytest.approx(11.5), "seconds")

    sl1 = samplelogs_instances[1]
    sl1.insert.assert_any_call("slice_start", pytest.approx(11.5), "seconds")
    sl1.insert.assert_any_call("slice_end", pytest.approx(21.5), "seconds")


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
