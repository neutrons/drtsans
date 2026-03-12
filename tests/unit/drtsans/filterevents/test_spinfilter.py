import pytest
from unittest.mock import MagicMock, patch

from drtsans.filterevents.spinfilter import SpinFilter, extract_times, create_table
from drtsans.polarization import (
    PV_POLARIZER,
    PV_ANALYZER,
    PV_POLARIZER_FLIPPER,
    PV_ANALYZER_FLIPPER,
    PV_POLARIZER_VETO,
    PV_ANALYZER_VETO,
)


def _make_workspace_group(xs_ids):
    """Create a mock workspace group whose items carry a cross_section_id log."""
    slices = [MagicMock(**{"getComment.return_value": xs_id}) for xs_id in xs_ids]
    group = MagicMock()
    group.getNumberOfEntries.return_value = len(slices)
    group.getItem.side_effect = lambda n: slices[n]
    return group


def _make_samplelogs(xs_id):
    """Return a SampleLogs mock that exposes cross_section_id."""
    sl = MagicMock()
    sl.__contains__ = MagicMock(side_effect=lambda key: key == "cross_section_id")
    sl.get = MagicMock(return_value=MagicMock(value=xs_id))
    return sl


def _make_device_sample_logs(has_polarizer=True, has_analyzer=True):
    """Return a SampleLogs mock that reports polarizer/analyzer device presence."""
    sl = MagicMock()
    sl.__contains__ = MagicMock(side_effect=lambda key: key in (PV_POLARIZER, PV_ANALYZER))
    sl.get = MagicMock(
        side_effect=lambda key, _: MagicMock(
            value=1 if ((key == PV_POLARIZER and has_polarizer) or (key == PV_ANALYZER and has_analyzer)) else 0
        )
    )
    return sl


def test_extract_times_multiple_time_values():
    """extract_times returns one tuple per timestamp."""
    times = [100, 200, 300]
    result = extract_times(times, device_on=True, is_polarizer=True)
    assert result == [
        (100, True, [True, False, False, False]),
        (200, True, [True, False, False, False]),
        (300, True, [True, False, False, False]),
    ]


def test_extract_times_different_device_states():
    """extract_times correctly sets all four device flags."""
    times = [100]
    device_on = True
    is_polarizer = True
    is_analyzer = True
    is_polarizer_veto = True
    is_analyzer_veto = True
    result = extract_times(times, device_on, is_polarizer, is_analyzer, is_polarizer_veto, is_analyzer_veto)
    assert result == [(100, True, [True, True, True, True])]


def test_extract_times_empty_list_of_times():
    """extract_times returns an empty list when given no timestamps."""
    times = []
    device_on = True
    result = extract_times(times, device_on)
    assert result == []


def test_create_table():
    """
    create_table discards the interval predating start_time and produces
    a first row starting at ~0 s for the Off_Off cross-section.
    """
    changes = [
        (1114103915483902567, False, [True, False, False, False]),
        (1114104949939068067, True, [False, False, True, False]),
        (1114104949966130867, True, [True, False, False, False]),
        (1114104949966130867, True, [True, False, False, False]),
        (1114104949968047467, False, [False, False, True, False]),
        (1114105024775176867, True, [False, False, True, False]),
        (1114105024785293967, False, [True, False, False, False]),
        (1114105024785293967, False, [True, False, False, False]),
        (1114105024785295267, False, [False, False, True, False]),
    ]
    table = create_table(changes, start_time=1114104876000000000, has_polarizer=True, has_analyzer=False)
    assert table.row(0) == {
        "start": pytest.approx(0.0, abs=0.1),
        "stop": pytest.approx(73.9, abs=0.1),
        "target": "Off_Off",
    }


@patch("drtsans.filterevents.spinfilter.workspace_handle")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_spin_filter_inject_metadata_common_fields(mock_samplelogs_cls, mock_workspace_handle):
    """inject_metadata inserts slice and number_of_slices for every cross-section."""
    xs_ids = ["On_On", "Off_Off"]
    mock_workspace_handle.return_value = _make_workspace_group(xs_ids)
    samplelogs_instances = [_make_samplelogs(xs_id) for xs_id in xs_ids]
    mock_samplelogs_cls.side_effect = samplelogs_instances

    SpinFilter(MagicMock()).inject_metadata("output_ws")

    sl0, sl1 = samplelogs_instances
    sl0.insert.assert_any_call("slice", 1)
    sl0.insert.assert_any_call("number_of_slices", 2)
    sl1.insert.assert_any_call("slice", 2)
    sl1.insert.assert_any_call("number_of_slices", 2)


@patch("drtsans.filterevents.spinfilter.workspace_handle")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_spin_filter_inject_metadata_polarization_fields(mock_samplelogs_cls, mock_workspace_handle):
    """inject_metadata inserts slice_parameter, cross_section, has_polarizer, has_analyzer."""
    xs_ids = ["On_On", "On_Off"]
    mock_workspace_handle.return_value = _make_workspace_group(xs_ids)
    samplelogs_instances = [_make_samplelogs(xs_id) for xs_id in xs_ids]
    mock_samplelogs_cls.side_effect = samplelogs_instances

    spin_filter = SpinFilter(MagicMock())
    spin_filter._has_polarizer = True
    spin_filter._has_analyzer = True
    spin_filter.inject_metadata("output_ws")

    sl0, sl1 = samplelogs_instances
    sl0.insert.assert_any_call("slice_parameter", "polarization_state")
    sl0.insert.assert_any_call("cross_section", "On_On")
    sl0.insert.assert_any_call("has_polarizer", 1)
    sl0.insert.assert_any_call("has_analyzer", 1)
    sl1.insert.assert_any_call("cross_section", "On_Off")


@patch("drtsans.filterevents.spinfilter.workspace_handle")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_spin_filter_inject_metadata_unknown_cross_section(mock_samplelogs_cls, mock_workspace_handle):
    """inject_metadata falls back to 'Unknown' when cross_section_id log is absent."""
    mock_workspace_handle.return_value = _make_workspace_group(["On_On"])
    sl = MagicMock()
    sl.__contains__ = MagicMock(return_value=False)  # cross_section_id not present
    mock_samplelogs_cls.return_value = sl

    SpinFilter(MagicMock()).inject_metadata("output_ws")

    sl.insert.assert_any_call("cross_section", "Unknown")


@patch("drtsans.filterevents.spinfilter.workspace_handle")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_spin_filter_inject_metadata_no_devices(mock_samplelogs_cls, mock_workspace_handle):
    """inject_metadata records has_polarizer=0, has_analyzer=0 when no devices are present."""
    mock_workspace_handle.return_value = _make_workspace_group(["Off_Off"])
    samplelogs_instances = [_make_samplelogs("Off_Off")]
    mock_samplelogs_cls.side_effect = samplelogs_instances

    spin_filter = SpinFilter(MagicMock())
    spin_filter._has_polarizer = False
    spin_filter._has_analyzer = False
    spin_filter.inject_metadata("output_ws")

    sl = samplelogs_instances[0]
    sl.insert.assert_any_call("has_polarizer", 0)
    sl.insert.assert_any_call("has_analyzer", 0)


# ---------------------------------------------------------------------------
# SpinFilter.__init__
# ---------------------------------------------------------------------------


def test_spin_filter_init_defaults():
    """__init__ stores PV names and initialises device-presence flags to False."""
    ws = MagicMock()
    sf = SpinFilter(ws)
    assert sf.workspace == str(ws)
    assert sf.pv_polarizer_state == PV_POLARIZER_FLIPPER
    assert sf.pv_analyzer_state == PV_ANALYZER_FLIPPER
    assert sf.pv_polarizer_veto == PV_POLARIZER_VETO
    assert sf.pv_analyzer_veto == PV_ANALYZER_VETO
    assert sf.check_devices is True
    assert sf._has_polarizer is False
    assert sf._has_analyzer is False


def test_spin_filter_init_custom_pv_names():
    """__init__ stores custom PV names when supplied."""
    sf = SpinFilter(
        MagicMock(),
        pv_polarizer_state="pol_state",
        pv_analyzer_state="ana_state",
        pv_polarizer_veto="pol_veto",
        pv_analyzer_veto="ana_veto",
        check_devices=False,
    )
    assert sf.pv_polarizer_state == "pol_state"
    assert sf.pv_analyzer_state == "ana_state"
    assert sf.pv_polarizer_veto == "pol_veto"
    assert sf.pv_analyzer_veto == "ana_veto"
    assert sf.check_devices is False


# ---------------------------------------------------------------------------
# SpinFilter.generate_filter
# ---------------------------------------------------------------------------


@patch("drtsans.filterevents.spinfilter.create_table")
@patch("drtsans.filterevents.spinfilter.workspace_handle")
@patch.object(SpinFilter, "_build_change_list")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_generate_filter_returns_empty_dict_when_devices_present(
    mock_samplelogs_cls, mock_build, mock_workspace_handle, mock_create_table
):
    """generate_filter returns {} and sets splitter_workspace when devices are found."""
    mock_samplelogs_cls.return_value = _make_device_sample_logs(has_polarizer=True, has_analyzer=True)
    mock_build.return_value = [(100, True, [True, False, False, False])]
    mock_workspace_handle.return_value.run.return_value.startTime.return_value.totalNanoseconds.return_value = 0
    mock_create_table.return_value = MagicMock()

    sf = SpinFilter(MagicMock())
    result = sf.generate_filter()

    assert result == {}
    assert sf._has_polarizer is True
    assert sf._has_analyzer is True
    mock_create_table.assert_called_once()


@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_generate_filter_returns_none_when_no_devices(mock_samplelogs_cls):
    """generate_filter returns None when neither polarizer nor analyzer is present."""
    mock_samplelogs_cls.return_value = _make_device_sample_logs(has_polarizer=False, has_analyzer=False)

    sf = SpinFilter(MagicMock())
    result = sf.generate_filter()

    assert result is None
    assert sf._has_polarizer is False
    assert sf._has_analyzer is False


@patch.object(SpinFilter, "_build_change_list")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_generate_filter_returns_none_when_empty_change_list(mock_samplelogs_cls, mock_build):
    """generate_filter returns None when the change list is empty."""
    mock_samplelogs_cls.return_value = _make_device_sample_logs(has_polarizer=True, has_analyzer=False)
    mock_build.return_value = []

    sf = SpinFilter(MagicMock())
    result = sf.generate_filter()

    assert result is None


@patch.object(SpinFilter, "_build_change_list")
@patch("drtsans.filterevents.spinfilter.SampleLogs")
def test_generate_filter_check_devices_false_assumes_both_present(mock_samplelogs_cls, mock_build):
    """When check_devices=False both polarizer and analyzer are assumed present."""
    mock_samplelogs_cls.return_value = _make_device_sample_logs(has_polarizer=False, has_analyzer=False)
    mock_build.return_value = []  # short-circuit before create_table

    sf = SpinFilter(MagicMock(), check_devices=False)
    sf.generate_filter()

    assert sf._has_polarizer is True
    assert sf._has_analyzer is True


# ---------------------------------------------------------------------------
# SpinFilter._build_change_list
# ---------------------------------------------------------------------------


@patch.object(SpinFilter, "_extract_veto_changes")
@patch.object(SpinFilter, "_extract_device_changes")
def test_build_change_list_only_polarizer(mock_device, mock_veto):
    """_build_change_list queries polarizer logs only when analyzer is absent."""
    mock_device.return_value = [(200, True, [True, False, False, False])]
    mock_veto.return_value = [(100, True, [False, False, True, False])]

    sf = SpinFilter(MagicMock())
    sf._has_polarizer = True
    sf._has_analyzer = False
    result = sf._build_change_list()

    mock_device.assert_called_once_with(sf.pv_polarizer_state, is_polarizer=True)
    mock_veto.assert_called_once_with(sf.pv_polarizer_veto, is_polarizer_veto=True)
    assert result == sorted(result, key=lambda x: x[0])


@patch.object(SpinFilter, "_extract_veto_changes")
@patch.object(SpinFilter, "_extract_device_changes")
def test_build_change_list_both_devices(mock_device, mock_veto):
    """_build_change_list queries both polarizer and analyzer device logs."""
    mock_device.return_value = [(300, True, [True, False, False, False])]
    mock_veto.return_value = []

    sf = SpinFilter(MagicMock())
    sf._has_polarizer = True
    sf._has_analyzer = True
    sf._build_change_list()

    assert mock_device.call_count == 2
    mock_device.assert_any_call(sf.pv_polarizer_state, is_polarizer=True)
    mock_device.assert_any_call(sf.pv_analyzer_state, is_analyzer=True)


@patch.object(SpinFilter, "_extract_veto_changes")
@patch.object(SpinFilter, "_extract_device_changes")
def test_build_change_list_result_is_sorted(mock_device, mock_veto):
    """_build_change_list returns changes sorted by timestamp."""
    mock_device.side_effect = [
        [(300, True, [True, False, False, False])],  # polarizer
        [(100, False, [False, True, False, False])],  # analyzer
    ]
    mock_veto.return_value = []

    sf = SpinFilter(MagicMock())
    sf._has_polarizer = True
    sf._has_analyzer = True
    result = sf._build_change_list()

    timestamps = [t for t, _, _ in result]
    assert timestamps == sorted(timestamps)


# ---------------------------------------------------------------------------
# SpinFilter._extract_device_changes
# ---------------------------------------------------------------------------


@patch("drtsans.filterevents.spinfilter.GenerateEventsFilter")
@patch("drtsans.filterevents.spinfilter.mtd")
def test_extract_device_changes_calls_generate_events_filter_twice(mock_mtd, mock_gef):
    """_extract_device_changes calls GenerateEventsFilter once for ON and once for OFF."""
    splitws = MagicMock()
    splitws.toDict.return_value = {"start": [100], "stop": [200]}
    mock_gef.return_value = (splitws, MagicMock())
    mock_mtd.unique_hidden_name.return_value = "_hidden"

    sf = SpinFilter(MagicMock())
    sf._extract_device_changes("some_log", is_polarizer=True)

    assert mock_gef.call_count == 2
    on_kwargs = mock_gef.call_args_list[0].kwargs
    assert on_kwargs["MinimumLogValue"] == pytest.approx(0.99)
    assert on_kwargs["MaximumLogValue"] == pytest.approx(1.01)
    off_kwargs = mock_gef.call_args_list[1].kwargs
    assert off_kwargs["MinimumLogValue"] == pytest.approx(-0.01)
    assert off_kwargs["MaximumLogValue"] == pytest.approx(0.01)


@patch("drtsans.filterevents.spinfilter.GenerateEventsFilter")
@patch("drtsans.filterevents.spinfilter.mtd")
def test_extract_device_changes_interprets_start_stop_correctly(mock_mtd, mock_gef):
    """
    ON filter: start→device_on=True, stop→device_on=False.
    OFF filter: start→device_on=False, stop→device_on=True.
    """
    on_ws = MagicMock()
    on_ws.toDict.return_value = {"start": [100], "stop": [200]}
    off_ws = MagicMock()
    off_ws.toDict.return_value = {"start": [300], "stop": [400]}
    mock_gef.side_effect = [(on_ws, MagicMock()), (off_ws, MagicMock())]
    mock_mtd.unique_hidden_name.return_value = "_hidden"

    sf = SpinFilter(MagicMock())
    result = sf._extract_device_changes("some_log", is_polarizer=True)

    assert (100, True, [True, False, False, False]) in result  # ON start  → on=True
    assert (200, False, [True, False, False, False]) in result  # ON stop   → on=False
    assert (300, False, [True, False, False, False]) in result  # OFF start → on=False
    assert (400, True, [True, False, False, False]) in result  # OFF stop  → on=True


# ---------------------------------------------------------------------------
# SpinFilter._extract_veto_changes
# ---------------------------------------------------------------------------


@patch("drtsans.filterevents.spinfilter.GenerateEventsFilter")
@patch("drtsans.filterevents.spinfilter.mtd")
def test_extract_veto_changes_calls_generate_events_filter_once(mock_mtd, mock_gef):
    """_extract_veto_changes calls GenerateEventsFilter exactly once (veto ON only)."""
    splitws = MagicMock()
    splitws.toDict.return_value = {"start": [500], "stop": [600]}
    mock_gef.return_value = (splitws, MagicMock())
    mock_mtd.unique_hidden_name.return_value = "_hidden"

    sf = SpinFilter(MagicMock())
    sf._extract_veto_changes("veto_log", is_polarizer_veto=True)

    assert mock_gef.call_count == 1
    call_kwargs = mock_gef.call_args.kwargs
    assert call_kwargs["MinimumLogValue"] == pytest.approx(0.99)
    assert call_kwargs["MaximumLogValue"] == pytest.approx(1.01)


@patch("drtsans.filterevents.spinfilter.GenerateEventsFilter")
@patch("drtsans.filterevents.spinfilter.mtd")
def test_extract_veto_changes_interprets_start_stop_correctly(mock_mtd, mock_gef):
    """Veto start→device_on=True (veto active), veto stop→device_on=False (lifted)."""
    splitws = MagicMock()
    splitws.toDict.return_value = {"start": [500], "stop": [600]}
    mock_gef.return_value = (splitws, MagicMock())
    mock_mtd.unique_hidden_name.return_value = "_hidden"

    sf = SpinFilter(MagicMock())
    result = sf._extract_veto_changes("veto_log", is_polarizer_veto=True)

    assert (500, True, [False, False, True, False]) in result
    assert (600, False, [False, False, True, False]) in result


# ---------------------------------------------------------------------------
# SpinFilter.apply_filter
# ---------------------------------------------------------------------------


@patch("drtsans.filterevents.spinfilter.GroupWorkspaces")
@patch.object(SpinFilter, "generate_filter")
def test_apply_filter_groups_raw_workspace_when_no_devices(mock_generate, mock_group):
    """apply_filter groups the raw workspace when generate_filter returns None."""
    mock_generate.return_value = None

    sf = SpinFilter("raw_ws")
    sf.apply_filter("output_ws")

    mock_group.assert_called_once_with(["raw_ws"], OutputWorkspace="output_ws")


@patch("drtsans.filterevents.spinfilter.GroupWorkspaces")
@patch.object(SpinFilter, "generate_filter")
def test_apply_filter_groups_raw_workspace_when_empty_splitter(mock_generate, mock_group):
    """apply_filter groups the raw workspace when the splitter table is empty."""
    mock_generate.return_value = {}
    splitter = MagicMock()
    splitter.rowCount.return_value = 0

    sf = SpinFilter("raw_ws")
    sf.splitter_workspace = splitter
    sf.apply_filter("output_ws")

    mock_group.assert_called_once_with(["raw_ws"], OutputWorkspace="output_ws")


@patch("drtsans.filterevents.spinfilter.AnalysisDataService")
@patch("drtsans.filterevents.spinfilter.AddSampleLog")
@patch("drtsans.filterevents.spinfilter.FilterEvents")
@patch("drtsans.filterevents.spinfilter.mtd")
@patch.object(SpinFilter, "generate_filter")
def test_apply_filter_calls_filter_events_and_adds_sample_log(
    mock_generate, mock_mtd, mock_filter_events, mock_add_log, mock_ads
):
    """apply_filter calls FilterEvents and adds cross_section_id to each output workspace."""
    mock_generate.return_value = {}
    splitter = MagicMock()
    splitter.rowCount.return_value = 2

    ws_on = MagicMock()
    ws_on.__str__ = lambda _: "output_ws_On_On"
    ws_off = MagicMock()
    ws_off.__str__ = lambda _: "output_ws_Off_Off"
    group = MagicMock()
    group.__iter__ = MagicMock(return_value=iter([ws_on, ws_off]))
    mock_filter_events.return_value = [MagicMock(), MagicMock(), group]
    mock_mtd.unique_hidden_name.return_value = "_corr"

    sf = SpinFilter("raw_ws")
    sf.splitter_workspace = splitter
    sf.apply_filter("output_ws")

    mock_filter_events.assert_called_once()
    mock_ads.remove.assert_called_once_with("_corr")
    mock_add_log.assert_any_call(Workspace=ws_on, LogName="cross_section_id", LogText="On_On")
    mock_add_log.assert_any_call(Workspace=ws_off, LogName="cross_section_id", LogText="Off_Off")
    ws_on.setComment.assert_called_once_with("On_On")
    ws_on.setTitle.assert_called_once_with("On_On")
    ws_off.setComment.assert_called_once_with("Off_Off")
    ws_off.setTitle.assert_called_once_with("Off_Off")
