# standard imports
import h5py
import os
import pytest
import tempfile

# third party imports
from mantid.kernel import Int64TimeSeriesProperty
from mantid.simpleapi import CreateWorkspace

# package imports
from drtsans.load import _insert_periodic_timeslice_log, __monitor_counts, resolve_slicing
from drtsans.samplelogs import SampleLogs


def test_monitor_counts():
    # Create HDF5 file with empty 'entry' group
    _, filename = tempfile.mkstemp(suffix=".h5")
    f = h5py.File(filename, "w")
    f.create_group("entry")
    f.close()
    # Assert we cannot read monitor counts
    with pytest.raises(RuntimeError) as except_info:
        __monitor_counts(filename)
    assert "does not contain /entry/" in str(except_info.value)
    # Append a monitor entry to the file
    f = h5py.File(filename, "a")
    group = f["entry"].create_group("monitor1")
    data_set = group.create_dataset("total_counts", (1,), dtype="i")
    total_counts = 42
    data_set[0] = total_counts
    f.close()
    # Assert the monitor counts
    assert __monitor_counts(filename) == total_counts
    os.remove(filename)


def test_periodic_timeslice_log(temp_workspace_name):
    # prepare the input workspace
    workspace = CreateWorkspace(dataX=[1.0], dataY=[42.0], OutputWorkspace=temp_workspace_name())
    sample_logs = SampleLogs(workspace)
    sample_logs.insert("duration", value=3600.0, unit="second")
    sample_logs.insert("run_start", value="2000-01-01T00:00:00", unit="")

    # insert the periodic log starting 42 seconds after run_start and having values from zero up to 60/1.5 -1 = 39
    _insert_periodic_timeslice_log(
        workspace, name="periodic_log", time_interval=1.5, time_period=60.0, time_offset=42.0
    )

    sample_logs = SampleLogs(workspace)
    assert "periodic_log" in sample_logs
    log = sample_logs["periodic_log"]
    assert isinstance(log, Int64TimeSeriesProperty)
    assert log.size() == (3600.0 - 42) / 1.5
    assert log.firstTime().toISO8601String() == "2000-01-01T00:00:42"
    assert "2000-01-01T00:59:58.5" in log.lastTime().toISO8601String()
    assert log.firstValue() == 0
    # the last period starts 42 seconds before the end of the run, thus (60 - 42) / 1.5 -1 = 11
    assert log.lastValue() == 11
    assert max(log.value) == 39


def test_resolve_slicing():
    options = {"useTimeSlice": True, "useLogSlice": False, "sample": {"runNumber": "12345"}}
    assert resolve_slicing(options) == (True, False)

    with pytest.raises(ValueError) as except_info:
        options["useLogSlice"] = True
        resolve_slicing(options)
    assert "Can't do both time and log slicing" in str(except_info.value)

    options = {"useTimeSlice": True, "useLogSlice": False, "sample": {"runNumber": "1,2,3"}}
    with pytest.raises(ValueError) as except_info:
        resolve_slicing(options)
    assert "Can't do slicing on summed data sets" in str(except_info.value)


if __name__ == "__main__":
    pytest.main([__file__])
