from unittest.mock import MagicMock, patch

from drtsans.filterevents.basefilter import FilterStrategy


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
