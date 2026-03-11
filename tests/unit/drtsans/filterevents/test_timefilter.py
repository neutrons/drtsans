import pytest
from unittest.mock import MagicMock, patch

from drtsans.filterevents.timefilter import TimeIntervalFilter, PeriodicTimeFilter


def _make_workspace_group(comments):
    # help function to create a mock workspace group with specified comments for each slice
    slices = [MagicMock(**{"getComment.return_value": c}) for c in comments]
    group = MagicMock()
    group.getNumberOfEntries.return_value = len(slices)
    group.getItem.side_effect = lambda n: slices[n]
    return group


def test_time_interval_filter_generate_filter():
    # default offset
    params = TimeIntervalFilter(MagicMock(), time_interval=10.0).generate_filter()
    assert params["TimeInterval"] == 10.0
    assert params["UnitOfTime"] == "Seconds"
    assert params["StartTime"] == "0.0"
    # explicit offset
    params = TimeIntervalFilter(MagicMock(), time_interval=5.0, time_offset=3.0).generate_filter()
    assert params["StartTime"] == "3.0"
    assert params["TimeInterval"] == 5.0


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.mtd")
@patch("drtsans.filterevents.timefilter.mtd")
def test_time_interval_inject_metadata(mock_timefilter_mtd, mock_base_mtd, mock_samplelogs_cls):
    run_start_ns = 1_000_000_000
    comments = ["slice_0", "slice_1"]
    mock_base_mtd.__getitem__.return_value = _make_workspace_group(comments)

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


# ------------------
# PeriodicTimeFilter
# ------------------


@patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log")
def test_periodic_time_filter_generate_filter(mock_create_log):
    params = PeriodicTimeFilter(MagicMock(), time_interval=5.0, time_period=60.0).generate_filter()
    assert params["LogName"] == "periodic_time_slicing"
    assert params["LogValueInterval"] == 1


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.mtd")
@patch("drtsans.filterevents.timefilter.PeriodicTimeFilter._create_periodic_log")
def test_periodic_time_inject_metadata(mock_create_log, mock_mtd, mock_samplelogs_cls):
    # Comment format parsed by the regex: "...From.X.To.Y.Value..."
    comments = [
        "Splitter.From.0.To.1.Value.0",
        "Splitter.From.1.To.2.Value.1",
    ]
    mock_mtd.__getitem__.return_value = _make_workspace_group(comments)
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
