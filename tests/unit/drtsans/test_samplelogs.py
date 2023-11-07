import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn
from mantid.simpleapi import LoadNexusProcessed
from mantid.api import Run

from drtsans.samplelogs import SampleLogs, periodic_index_log


def test_periodic_index_log():
    period = 1.0  # one second
    interval = 1.0 / 30  # 30 Hz
    duration = 3600.5  # one hour and half a second
    offset = 60.0  # one minute
    run_start = "2023-09-02T12:00:00"  # beginning of the run

    log = periodic_index_log(period, interval, duration, run_start, offset)
    assert log.name == "periodic_index"
    assert log.size() == (3600.5 - 60) * 30
    assert log.firstTime().toISO8601String() == "2023-09-02T12:01:00"
    assert log.firstValue() == 0
    assert "2023-09-02T13:00:00.466" in log.lastTime().toISO8601String()
    assert log.lastValue() == 14


class TestSampleLogs:
    @pytest.mark.datarepo
    def test_init(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        w = LoadNexusProcessed(test_file, OutputWorkspace=clean_workspace("test_init_w"))
        r = w.getRun()
        for other in [w, r]:
            sl = SampleLogs(other)
            assert isinstance(sl._run, Run)

    @pytest.mark.datarepo
    def test_getitem(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        ws = LoadNexusProcessed(Filename=test_file)
        clean_workspace(ws)
        sl = SampleLogs(ws)
        assert_almost_equal(sl["Phase1"].value.mean(), 22444, decimal=0)

        with pytest.raises(KeyError):
            sl["nonexistantlog"].value
            assert False, "Should have failed \"sl['nonexistantlog'].value\""

    @pytest.mark.datarepo
    def test_getattr(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        ws = LoadNexusProcessed(Filename=test_file)
        clean_workspace(ws)
        sl = SampleLogs(ws)
        assert_almost_equal(sl.Phase1.value.mean(), 22444, decimal=0)

        with pytest.raises(AttributeError):
            sl.nonexistantlog.value
            assert False, 'Should have failed "sl.nonexistantlog.value"'

    @pytest.mark.datarepo
    def test_insert(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        ws = LoadNexusProcessed(test_file)
        clean_workspace(ws)
        sl = SampleLogs(ws)
        sl.insert("string_log", "log value")
        assert sl.string_log.value, "log value"
        assert not sl.string_log.units

        units = "super awesome units"
        sl.insert("int_log", 42, units)
        assert sl.int_log.value == 42
        assert sl.int_log.units == units

        euler = 2.7182818284590452353602874713527
        units = "even more awesome units"
        sl.insert("float_log", euler, units)
        assert sl.float_log.units == units
        assert sl.float_log.value == euler

        values = list(range(1, 9))
        units = "most awesomest units ever"
        sl.insert("array_log", values, units)
        assert sl.array_log.units == units
        # this seems like the wrong value from mantid
        assert sl.array_log.value == values[0]

        period = 1.0  # one second
        interval = 1.0 / 30  # 30 Hz
        duration = 3600.5  # one hour and half a second
        offset = 60.0  # one minute
        run_start = "2023-09-02T12:00:00"  # beginning of the run
        log = periodic_index_log(period, interval, duration, run_start, offset, name="periodic")
        sl.insert(name=log.name, value=log)
        assert sl.periodic.name == "periodic"
        assert sl.periodic.size() == (3600.5 - 60) * 30
        assert sl.periodic.firstTime().toISO8601String() == "2023-09-02T12:01:00"
        assert sl.periodic.firstValue() == 0
        assert "2023-09-02T13:00:00.466" in sl.periodic.lastTime().toISO8601String()
        assert sl.periodic.lastValue() == 14


if __name__ == "__main__":
    pytest.main([__file__])
