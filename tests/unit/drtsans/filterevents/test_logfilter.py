from unittest.mock import MagicMock, patch
from drtsans.filterevents.logfilter import LogValueFilter


def _make_workspace_group(comments):
    # Create a mock workspace group with the specified comments for each slice
    slices = [MagicMock(**{"getComment.return_value": c}) for c in comments]
    group = MagicMock()
    group.getNumberOfEntries.return_value = len(slices)
    group.getItem.side_effect = lambda n: slices[n]
    return group


def _make_samplelogs(units):
    """Return a SampleLogs mock whose log lookup returns the given units string."""
    sl = MagicMock()
    sl.__contains__ = MagicMock(return_value=(units is not None))
    sl.__getitem__ = MagicMock(return_value=MagicMock(units=units or ""))
    return sl


def test_log_value_filter_generate_filter():
    params = LogValueFilter(MagicMock(), log_name="SampleTemp", log_value_interval=5.0).generate_filter()
    assert params == {"LogName": "SampleTemp", "LogValueInterval": 5.0}


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.mtd")
def test_log_value_inject_metadata_with_units(mock_mtd, mock_samplelogs_cls):
    comments = [
        "Splitter.From.10.To.15.Value.0",
        "Splitter.From.15.To.20.Value.1",
    ]
    mock_mtd.__getitem__.return_value = _make_workspace_group(comments)
    samplelogs_instances = [_make_samplelogs(units="K"), _make_samplelogs(units="K")]
    mock_samplelogs_cls.side_effect = samplelogs_instances
    LogValueFilter(MagicMock(), log_name="SampleTemp", log_value_interval=5.0).inject_metadata("output_ws")
    sl0, sl1 = samplelogs_instances
    sl0.insert.assert_any_call("slice_parameter", "SampleTemp")
    sl0.insert.assert_any_call("slice_interval", 5.0)
    sl0.insert.assert_any_call("slice_start", 10.0, "K")
    sl0.insert.assert_any_call("slice_end", 15.0, "K")
    sl1.insert.assert_any_call("slice_start", 15.0, "K")
    sl1.insert.assert_any_call("slice_end", 20.0, "K")


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.mtd")
def test_log_value_inject_metadata_without_units(mock_mtd, mock_samplelogs_cls):
    comments = ["Splitter.From.3.To.6.Value.0"]
    mock_mtd.__getitem__.return_value = _make_workspace_group(comments)
    samplelogs_instances = [_make_samplelogs("")]  # present but no units
    mock_samplelogs_cls.side_effect = samplelogs_instances
    LogValueFilter(MagicMock(), log_name="ProtonCharge", log_value_interval=3.0).inject_metadata("output_ws")
    sl = samplelogs_instances[0]
    sl.insert.assert_any_call("slice_start", 3.0)
    sl.insert.assert_any_call("slice_end", 6.0)


@patch("drtsans.filterevents.basefilter.SampleLogs")
@patch("drtsans.filterevents.basefilter.mtd")
def test_log_value_inject_metadata_common_fields(mock_mtd, mock_samplelogs_cls):
    comments = [
        "Splitter.From.0.To.5.Value.0",
        "Splitter.From.5.To.10.Value.1",
    ]
    mock_mtd.__getitem__.return_value = _make_workspace_group(comments)
    samplelogs_instances = [_make_samplelogs(None), _make_samplelogs(None)]
    mock_samplelogs_cls.side_effect = samplelogs_instances
    LogValueFilter(MagicMock(), log_name="SampleTemp", log_value_interval=5.0).inject_metadata("output_ws")
    sl0, sl1 = samplelogs_instances
    sl0.insert.assert_any_call("slice", 1)
    sl0.insert.assert_any_call("number_of_slices", 2)
    sl0.insert.assert_any_call("slice_info", comments[0])
    sl1.insert.assert_any_call("slice", 2)
    sl1.insert.assert_any_call("number_of_slices", 2)
    sl1.insert.assert_any_call("slice_info", comments[1])
