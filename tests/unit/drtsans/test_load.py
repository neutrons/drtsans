# standard imports
import h5py
import os
import pytest
import tempfile
from unittest import mock

# third party imports
from mantid.kernel import FloatTimeSeriesProperty
from mantid.simpleapi import CreateWorkspace, LoadEmptyInstrument, mtd
import numpy as np
from numpy.testing import assert_array_almost_equal

# package imports
from drtsans.load import _insert_periodic_timeslice_log, __monitor_counts, load_events, resolve_slicing
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


@pytest.mark.parametrize(
    "ID, time_interval",
    [
        ("integerNslices", 1.5),  # period is a multiple of time slice interval
        ("nonIntegerNslices", 8.0),  # ... or it is not
    ],
)
def test_periodic_timeslice_log(temp_workspace_name, ID, time_interval):
    # prepare the input workspace

    duration = 3600.0
    time_period = 60.0
    time_offset = 42.0

    workspace = CreateWorkspace(dataX=[1.0], dataY=[42.0], OutputWorkspace=temp_workspace_name())
    sample_logs = SampleLogs(workspace)
    sample_logs.insert("duration", value=duration, unit="second")
    sample_logs.insert("run_start", value="2000-01-01T00:00:00", unit="")

    # insert the periodic log starting 42 seconds after run_start
    # and having values from zero up to period / time_interval - 1
    _insert_periodic_timeslice_log(
        workspace, name="periodic_log", time_interval=time_interval, time_period=time_period, time_offset=time_offset
    )

    sample_logs = SampleLogs(workspace)
    log = sample_logs["periodic_log"]

    # common asserts
    assert "periodic_log" in sample_logs
    assert isinstance(log, FloatTimeSeriesProperty)
    assert log.firstTime().toISO8601String() == "2000-01-01T00:00:42"
    assert log.firstValue() == 0

    # case-wise asserts
    if ID == "integerNslices":
        # there are (duration - offset) / timesliceInterval entries
        assert log.size() == (duration - time_offset) / time_interval
        assert "2000-01-01T00:59:58.5" in log.lastTime().toISO8601String()
        # the last period starts 42 seconds before the end of the run, thus
        assert log.lastValue() == (time_period - time_offset) / time_interval - 1
        # indices range [0 : 60 / 1.5 = 40 : 1]
        assert max(log.value) == (time_period / time_interval) - 1

    elif ID == "nonIntegerNslices":
        n_periods = np.floor((duration - time_offset) / time_period)
        n_intervals_period = np.ceil(time_period / time_interval)
        n_intervals_remainder = np.ceil(((duration - time_offset) % time_period) / time_interval)

        assert log.size() == n_periods * n_intervals_period + n_intervals_remainder
        assert "2000-01-01T00:59:58" in log.lastTime().toISO8601String()
        # and again, [0 : 60 / 8 = 7.5 : 1]
        assert log.lastValue() == np.floor((60 - 42) / 8)
        assert max(log.value) == np.floor(60 / 8)


def test_resolve_slicing():
    options = {"configuration": {"useTimeSlice": True, "useLogSlice": False}, "sample": {"runNumber": "12345"}}
    assert resolve_slicing(options) == (True, False)

    with pytest.raises(ValueError) as except_info:
        options["configuration"]["useLogSlice"] = True
        resolve_slicing(options)
    assert "Can't do both time and log slicing" in str(except_info.value)

    options = {"configuration": {"useTimeSlice": True, "useLogSlice": False}, "sample": {"runNumber": "1,2,3"}}
    with pytest.raises(ValueError) as except_info:
        resolve_slicing(options)
    assert "Can't do slicing on summed data sets" in str(except_info.value)


class TestLoadEvents:
    def test_scale_components(self, temp_workspace_name, tmp_path):
        r"""Test application of Mantid algorithm ScaleInstrumentComponent inside src.drtsans.load.load_events

        Create an empty instrument workspace for each instrument, then mock all the necessary
        functions used in load_events() so that steps carrying out the loading of the event data
        ends up returning the empty instrument workspace.
        After that, apply ScaleInstrumentComponent to the workspace.
        After that, mock all other functions used in load_events(). The effect of these functions is not
        relevant to this test.
        """
        workspace = temp_workspace_name()
        scalings = {"detector1": [1.1, 1.2, 1.3]}

        # mock various functions used in load_events that are not relevant to this test
        with (
            mock.patch("drtsans.load.LoadEventNexus") as mock_load_event_nexus,
            mock.patch("drtsans.load.LoadEventAsWorkspace2D") as mock_load_event_as_workspace2d,
            mock.patch("drtsans.load.__monitor_counts") as mock_monitor_counts,
            mock.patch("drtsans.geometry.sample_detector_distance") as mock_sample_detector_distance,
        ):
            mock_monitor_counts.return_value = 42
            mock_sample_detector_distance.return_value = 0.42

            for instrument in ["CG2", "CG3", "EQSANS"]:
                # create empty file `runfile` only to fool the functionality inside load_events() that
                # verifies the existence of the input event data file
                runfile = tmp_path / f"{instrument}_12345.nxs.h5"
                runfile.touch()
                LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=workspace)
                mock_load_event_nexus.return_value = None
                mock_load_event_as_workspace2d.return_value = None
                load_events(run=str(runfile), scale_components=scalings, output_workspace=workspace)

                # check that the scalings [1.1, 1.2, 1.3] have been applied to panel 'detector1'
                compInfo = mtd[workspace].componentInfo()
                detectors = compInfo.detectorsInSubtree(compInfo.indexOfAny("detector1"))
                step_size = int(len(detectors) / 10)  # check 10 detectors only
                for i in range(0, len(detectors), step_size):
                    detector_id = detectors[i]
                    assert_array_almost_equal(compInfo.scaleFactor(int(detector_id)), scalings["detector1"])


if __name__ == "__main__":
    pytest.main([__file__])
