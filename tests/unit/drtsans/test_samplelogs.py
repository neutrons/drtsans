from os.path import join as pjn

import pytest
from mantid.api import Run
from mantid.simpleapi import LoadNexusProcessed
from numpy.testing import assert_almost_equal

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


@pytest.mark.parametrize(
    "interval, period, duration, offset, step, expected",
    [
        # integer number of intervals in period
        (1, 2, 10, 0, 1, 10),
        # non-integer number of intervals in period
        (0.01, 0.9901076314875992, 8371.0458984375, 0.6694618359281799, 1.0, 845401),
        (1.0, 9.0914642068604, 382.9853210449219, 1.1055757465580345, 1.0, 421),
        (1.0, 9.091311484270927, 380.5020751953125, 4.384090426034748, 1.0, 414),
        (1.0, 90.87363154577008, 227.791259765625, 63.92350945782564, 1.0, 164),
        (0.01, 0.990110528971648, 8353.4296875, 0.4783230672977612, 1.0, 843638),
        (0.1, 0.990107380210766, 7111.37744140625, 0.8365827919520669, 1.0, 71816),
        (0.1, 0.9901, 763.654052734375, 0.578, 1.0, 7707),
        (0.01, 0.09990091554886044, 228.7245635986328, 0.04903141722988207, 1.0, 22891),
        (0.01, 0.09990091554886044, 1903.1103515625, 0.04903141722988207, 1.0, 190495),
        (0.01, 0.09990100006257276, 7.999692916870117, 0.06871597503197445, 1.0, 794),
        (0.009, 0.09990099085909061, 8.049692153930664, 0.011064958593924914, 1.0, 966),
        (0.009, 0.09990099496711198, 38.19853591918945, 0.06973558286231382, 1.0, 4580),
        # interval = period | duration reads:
        # interval equals period, period factors duration
        #
        # interval !| period !| duration reads:
        # interval does not factor period, period does not factor duration
        # interval = period = duration
        (60, 60, 60, 0, 1, 1),
        # interval | period = duration
        (30, 60, 60, 0, 1, 2),
        # interval !| period = duration
        (40, 60, 60, 0, 1, 2),
        # interval = period | duration
        (30, 30, 60, 0, 1, 2),
        # interval | period | duration
        (15, 30, 60, 0, 1, 4),
        # interval | period !| duration
        (25, 50, 60, 0, 1, 3),
        # interval !| period | duration
        (25, 30, 60, 0, 1, 4),
        # interval !| period !| duration
        (40, 50, 60, 0, 1, 3),
        # interval | duration, (but does not factor period, covered above)
        (30, 40, 60, 0, 1, 3),
        # offset tests
        (60, 60, 60, 1, 1, 1),
        # edge case where entries must overtile, or there won't be enough to match times
        (0.31416253927, 3.1416253927, 11063.2470703125, 1.19011221899, 1, 38639),
    ],
    ids=[
        # integer number of intervals in period
        "period_count-exact-integer",
        # non-integer number of intervals in period
        "reduce_slice_157044",
        "reduce_slice_157054",
        "reduce_slice_157065",
        "reduce_slice_157076",
        "reduce_slice_157045",
        "reduce_slice_157046",
        "reduce_slice_157053",
        "reduce_slice_157106",
        "reduce_slice_157107",
        "reduce_slice_157110",
        "reduce_slice_157111",
        "reduce_slice_157112",
        "i=p|d",
        "i|p=d",
        "i!|p=d",
        "i=p|d",
        "i|p|d",
        "i|p!|d",
        "i!|p|d",
        "i!|p!|d",
        "i|d",
        "offset-duration",
        "require-overtiling",
    ],
)
def test_periodic_index_log_cases(interval, period, duration, offset, step, expected):
    run_start = "2000-01-01T00:00:00"
    if expected == "fail":
        with pytest.raises(AssertionError, match=r"(AssertionError: times and entries must have the same length)*"):
            result = periodic_index_log(period, interval, duration, run_start, offset, step)
    else:
        result = periodic_index_log(period, interval, duration, run_start, offset, step)
        assert result.size() == expected, f"Expected TimeSeriesProperty with size {expected}, got {result.size()}"


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

    @pytest.mark.datarepo
    def test_single_value(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        ws = LoadNexusProcessed(Filename=test_file)
        clean_workspace(ws)
        sl = SampleLogs(ws)
        assert isinstance(sl.single_value("run_start"), str)
        assert sl.single_value("sample_detector_distance") == sl.sample_detector_distance.value
        assert_almost_equal(sl.single_value("proton_charge"), sl._run.getStatistics("proton_charge").mean, decimal=0)
        with pytest.raises(RuntimeError):
            sl.single_value("log_does_not_exist")

    @pytest.mark.datarepo
    def test_records_json(self, datarepo_dir, clean_workspace):
        test_file = pjn(datarepo_dir.sans, "test_samplelogs", "EQSANS_92353_no_events.nxs")
        ws = LoadNexusProcessed(Filename=test_file)
        clean_workspace(ws)
        sl = SampleLogs(ws)
        records = sl.records_json(log="run_start", logs=["proton_charge", "log_does_not_exist"], default=None)
        assert records["run_start"] == "2017-12-19T21:38:25.356596666"
        assert_almost_equal(records["proton_charge"], 10064207.3, decimal=0)
        assert records["log_does_not_exist"] is None


if __name__ == "__main__":
    pytest.main([__file__])
