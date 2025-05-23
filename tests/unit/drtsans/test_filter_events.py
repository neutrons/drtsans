import pytest
from drtsans.filter_events import create_table, extract_times


def test_extract_times_multiple_time_values():
    """Test extract_times with a list of times"""
    times = [100, 200, 300]
    result = extract_times(times, is_start=True, is_polarizer=True)
    assert result == [
        (100, True, [True, False, False, False]),
        (200, True, [True, False, False, False]),
        (300, True, [True, False, False, False]),
    ]


def test_extract_times_different_device_states():
    """Test extract_times with different device states"""
    times = [100]
    is_start = True
    is_polarizer = True
    is_analyzer = True
    is_polarizer_veto = True
    is_analyzer_veto = True
    result = extract_times(times, is_start, is_polarizer, is_analyzer, is_polarizer_veto, is_analyzer_veto)
    assert result == [(100, True, [True, True, True, True])]


def test_extract_times_empty_list_of_times():
    """Test extract_times with an empty list of times"""
    times = []
    is_start = True
    result = extract_times(times, is_start)
    assert result == []


def test_create_table():
    """
    Test the creation of a table with input of a list containing an item (the first item)
    with a time stamp that predates the start time.
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
