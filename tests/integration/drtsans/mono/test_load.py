# Test load GPSANS and BIOSANS data
import pytest
import os
import numpy as np
from drtsans.filterevents.timefilter import TimeIntervalFilter
from drtsans.mono.load import load_and_split, load_events
from drtsans.mono.meta_data import get_sample_detector_offset
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import sample_detector_distance
from drtsans.load import move_instrument
from mantid.simpleapi import AddSampleLogMultiple, LoadEventNexus, mtd


def test_load_gpsans(datarepo_dir, clean_workspace):
    """Test load GPSANS data

    Returns
    -------

    """

    nexus_file_name = os.path.join(datarepo_dir.gpsans, "CG2_9165.nxs.h5")

    # Load data
    ws = load_events(
        nexus_file_name,
        output_workspace="gptest01",
        overwrite_instrument=True,
        detector_offset=0,
        sample_offset=0,
    )
    clean_workspace(ws)  # delete on test completion

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print("[TEST INFO] SampleToSi = {} mm".format(logs.find_log_with_units("CG2:CS:SampleToSi", unit="mm")))
    raw_sample_det_distance = sample_detector_distance(ws, unit="m", search_logs=False)
    print(
        "[TEST INFO] Sample to detector distance = {} /{} meter".format(
            raw_sample_det_distance,
            sample_detector_distance(ws, unit="m", log_key="sample_detector_distance", search_logs=True),
        )
    )

    # sample and detector offsets can only be retrieved from a loaded workspace
    # This is a technical debt
    sample_offset, detector_offset = get_sample_detector_offset(ws, "CG2:CS:SampleToSi", 0.0)

    assert sample_offset == -0.083
    assert detector_offset == -0.083

    # Move instrument
    # Move sample and detector
    ws = move_instrument(
        ws,
        sample_offset,
        detector_offset,
        is_mono=True,
        sample_si_name="CG2:CS:SampleToSi",
        si_window_to_nominal_distance=0,
    )

    # Verify
    new_sample_det_distance = sample_detector_distance(ws, unit="m", search_logs=False)
    print("[TEST INFO] Sample detector distance after moving = {} meter".format(new_sample_det_distance))
    print("[TEST INFO] Sample position = {}".format(ws.getInstrument().getSample().getPos()))

    assert new_sample_det_distance == raw_sample_det_distance


def test_load_biosans(datarepo_dir, clean_workspace):
    """Test load BIOSANS data

    Returns
    -------

    """

    nexus_file_name = os.path.join(datarepo_dir.biosans, "CG3_5705.nxs.h5")

    # Load data
    ws = load_events(
        nexus_file_name,
        output_workspace="biotest01",
        overwrite_instrument=True,
        detector_offset=0,
        sample_offset=0,
    )
    clean_workspace(ws)

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print("[TEST INFO] (Raw) sampleToSi = {} mm".format(logs.find_log_with_units("CG3:CS:SampleToSi", unit="mm")))
    raw_sample_det_distance = sample_detector_distance(ws)
    print(
        "[TEST INFO] (Raw) sample to detector distance = {} /{} meter".format(
            raw_sample_det_distance,
            sample_detector_distance(ws, log_key="sample_detector_distance", search_logs=True),
        )
    )

    # Calculate offset without any overwriting
    # sample and detector offsets can only be retrieved from a loaded workspace
    # This is a technical debt
    sample_offset, detector_offset = get_sample_detector_offset(ws, "CG3:CS:SampleToSi", 71.0 * 1e-3)
    print("[TEST INFO] Sample offset = {}, Detector offset = {}".format(sample_offset, detector_offset))

    # Verify: No sample offset from nominal position (origin)
    assert sample_offset == pytest.approx(0.0, 1e-12)
    # Verify: No sample offset
    assert detector_offset == pytest.approx(0.0, 1e-12)
    # Verify: sample position at (0., 0., 0.)
    sample_pos = np.array(ws.getInstrument().getSample().getPos())
    expected_sample_pos = np.array([0.0, 0.0, 0.0])
    np.testing.assert_allclose(sample_pos, expected_sample_pos, atol=1e-12)
    # Verify: sample detector distance is equal to 7.00000019 meter
    sample_det_distance_cal = sample_detector_distance(ws, unit="m", search_logs=False)
    sample_det_distance_meta = sample_detector_distance(ws, unit="mm", search_logs=True)
    assert sample_det_distance_cal == pytest.approx(sample_det_distance_meta * 1e-3, 1e-7)
    assert sample_det_distance_cal == pytest.approx(7.00000019, 1e-7)


def test_load_biosans_sample_off_nominal(datarepo_dir, clean_workspace):
    """Test load BIOSANS data with sample position off nominal position

    Returns
    -------

    """

    nexus_file_name = os.path.join(datarepo_dir.biosans, "CG3_5705.nxs.h5")

    # Load data
    ws = load_events(
        nexus_file_name,
        output_workspace="biotest01",
        overwrite_instrument=True,
        detector_offset=0,
        sample_offset=0,
    )
    clean_workspace(ws)

    # Verify: sample position at (0., 0., 0.)
    sample_pos = np.array(ws.getInstrument().getSample().getPos())
    expected_sample_pos = np.array([0.0, 0.0, 0.0])
    np.testing.assert_allclose(sample_pos, expected_sample_pos, atol=1e-12)
    # Verify: sample detector distance is equal to 7.00000019 meter
    sample_det_distance_cal = sample_detector_distance(ws, unit="m", search_logs=False)
    sample_det_distance_meta = sample_detector_distance(ws, unit="mm", search_logs=True)
    assert sample_det_distance_cal == pytest.approx(sample_det_distance_meta * 1e-3, 1e-7)
    assert sample_det_distance_cal == pytest.approx(7.00000019, 1e-7)

    # Second test on SampleToSi distance other than 71.00 mm
    # Simulate sample log CG3:CS:SampleToSi
    # Note: this is not overwriting SampleSiDistance
    test_sample_si_distance = 74.21
    AddSampleLogMultiple(
        Workspace=ws,
        LogNames="{}".format("CG3:CS:SampleToSi"),
        LogValues="{}".format(test_sample_si_distance),
        LogUnits="mm",
    )

    # Calculate offset without any overwriting
    sample_offset, detector_offset = get_sample_detector_offset(ws, "CG3:CS:SampleToSi", 71.0 * 1e-3)
    print("[TEST INFO 2] Sample offset = {}, Detector offset = {}".format(sample_offset, detector_offset))

    # Both sample and detector shall move toward souce (-Y direction) with (74.21 - 71.) = 3.21 mm
    assert sample_offset == pytest.approx(-0.00321, 1e-12)
    assert detector_offset == pytest.approx(-0.00321, 1e-12)

    # Move sample and detector
    ws = move_instrument(
        ws,
        sample_offset,
        detector_offset,
        is_mono=True,
        sample_si_name="CG3:CS:SampleToSi",
        si_window_to_nominal_distance=0,
    )

    # Verify: sample position at (0., 0., -0.00321)
    sample_pos = np.array(ws.getInstrument().getSample().getPos())
    expected_sample_pos = np.array([0.0, 0.0, -0.00321])
    np.testing.assert_allclose(sample_pos, expected_sample_pos, atol=1e-12)

    # Verify the sample detector distance which shall be same as raw meta data
    sample_det_distance_cal = sample_detector_distance(ws, unit="m", search_logs=False)
    sample_det_distance_meta = sample_detector_distance(ws, unit="mm", search_logs=True)
    assert sample_det_distance_cal == pytest.approx(sample_det_distance_meta * 1e-3, 1e-7)
    assert sample_det_distance_cal == pytest.approx(7.00000019, 1e-7)


def test_load_biosans_overwrite_swd(datarepo_dir, clean_workspace):
    """Test load BIOSANS data with overwriting sample Si window distance

    Returns
    -------

    """

    nexus_file_name = os.path.join(datarepo_dir.biosans, "CG3_5705.nxs.h5")

    # Load data
    ws = load_events(
        nexus_file_name,
        output_workspace="biotest02",
        overwrite_instrument=True,
        detector_offset=0,
        sample_offset=0,
    )
    clean_workspace(ws)

    # Calculate offset with overwriting to sample-detector-distance
    sample_offset, detector_offset = get_sample_detector_offset(
        ws, "CG3:CS:SampleToSi", 71.0 * 1e-3, overwrite_sample_si_distance=0.07421
    )
    print("[TEST INFO] Sample offset = {}, Detector offset = {}".format(sample_offset, detector_offset))

    # Move sample and detector
    ws = move_instrument(
        ws,
        sample_offset,
        detector_offset,
        is_mono=True,
        sample_si_name="CG3:CS:SampleToSi",
        si_window_to_nominal_distance=0.071,
    )

    # Verify: sample position at (0., 0., -0.00321) because SampleToSi is overwritten to 74.21 mm
    sample_pos = np.array(ws.getInstrument().getSample().getPos())
    expected_sample_pos = np.array([0.0, 0.0, -0.00321])
    np.testing.assert_allclose(sample_pos, expected_sample_pos, atol=1e-12)

    # Verify the sample detector distance shall be increased by 3.21mm due to the shift of sample position
    sample_det_distance_cal = sample_detector_distance(ws, unit="m", search_logs=False)
    assert sample_det_distance_cal == pytest.approx(7.00321, 1e-7)
    # verify the values from calculated and from meta data are identical
    sample_det_distance_meta = sample_detector_distance(ws, unit="mm", search_logs=True)
    assert sample_det_distance_cal == pytest.approx(sample_det_distance_meta * 1e-3, 1e-7)
    # verify that SampleToSi is overwritten to 74.21 mm
    logs = SampleLogs(ws)
    swd = logs.find_log_with_units("CG3:CS:SampleToSi", unit="mm")
    assert swd == pytest.approx(74.21, 1e-10)


def test_load_biosans_overwrite_sdd(datarepo_dir, clean_workspace):
    """Test load BIOSANS data with overwriting sample detector distance related meta data

    Returns
    -------

    """

    nexus_file_name = os.path.join(datarepo_dir.biosans, "CG3_5705.nxs.h5")

    # Load data
    ws = load_events(
        nexus_file_name,
        output_workspace="biotest02",
        overwrite_instrument=True,
        detector_offset=0,
        sample_offset=0,
    )
    clean_workspace(ws)

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print("[TEST INFO] SampleToSi = {} mm".format(logs.find_log_with_units("CG3:CS:SampleToSi", unit="mm")))
    raw_sample_det_distance = sample_detector_distance(ws)
    print(
        "[TEST INFO] Sample to detector distance = {} /{} meter".format(
            raw_sample_det_distance,
            sample_detector_distance(ws, log_key="sample_detector_distance", search_logs=True),
        )
    )

    # Calculate offset with overwriting to sample-detector-distance
    sample_offset, detector_offset = get_sample_detector_offset(
        ws, "CG3:CS:SampleToSi", 71.0 * 1e-3, overwrite_sample_detector_distance=7.1234
    )
    print("[TEST INFO] Sample offset = {}, Detector offset = {}".format(sample_offset, detector_offset))

    # Move sample and detector
    ws = move_instrument(
        ws,
        sample_offset,
        detector_offset,
        is_mono=True,
        sample_si_name="CG3:CS:SampleToSi",
        si_window_to_nominal_distance=0,
    )

    # Verify: sample position at (0., 0., 0.) because SampleToSi == 71 mm and not overwritten
    sample_pos = np.array(ws.getInstrument().getSample().getPos())
    expected_sample_pos = np.array([0.0, 0.0, 0.0])
    np.testing.assert_allclose(sample_pos, expected_sample_pos, atol=1e-12)

    # Verify the sample detector distance which shall be same as raw meta data
    sample_det_distance_cal = sample_detector_distance(ws, unit="m", search_logs=False)
    assert sample_det_distance_cal == pytest.approx(7.1234, 1e-7)
    # verify the values from calculated and from meta data are identical
    sample_det_distance_meta = sample_detector_distance(ws, unit="mm", search_logs=True)
    assert sample_det_distance_cal == pytest.approx(sample_det_distance_meta * 1e-3, 1e-7)


@pytest.mark.datarepo
def test_time_slice_metadata_with_load_time_cutoffs(datarepo_dir):
    """Slices emptied by FilterByTimeStart/FilterByTimeStop are deleted; the remaining slices keep their windows."""
    # Arrange: 300 s GPSANS run. Slices every 15 s from 5 s; only events between 40 s and 100 s are loaded,
    # so the slices 5-20 s and 20-35 s are empty and removed, as are all slices after 110 s.
    filename = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
    existing_workspaces = set(mtd.getObjectNames())

    try:
        # Act
        split_group = load_and_split(
            filename,
            sample_to_si_name="CG2:CS:SampleToSi",
            si_nominal_distance=0.0,
            time_interval=15.0,
            time_offset=5.0,
            output_workspace="time_slice_metadata_cutoffs",
            FilterByTimeStart=40.0,
            FilterByTimeStop=100.0,
        )

        # Assert
        assert split_group.getNumberOfEntries() == 5
        # the removed empty slices and their monitor slices are deleted, not only taken out of their groups
        monitor_group = mtd["time_slice_metadata_cutoffs_monitors"]
        kept = {split_group.name(), monitor_group.name()}
        kept |= set(split_group.getNames()) | set(monitor_group.getNames())
        created = set(mtd.getObjectNames()) - existing_workspaces
        assert {name for name in created if name.startswith("time_slice_metadata_cutoffs")} == kept
        for n, expected_start in enumerate([35.0, 50.0, 65.0, 80.0, 95.0]):
            slice_workspace = split_group.getItem(n)
            sample_logs = SampleLogs(slice_workspace)
            assert sample_logs.slice_start.value == pytest.approx(expected_start)
            assert sample_logs.slice_end.value == pytest.approx(expected_start + 15.0)
            # the first pulse of the slice lies inside its recorded window
            run_start = slice_workspace.run().startTime()
            first_pulse = (sample_logs["proton_charge"].firstTime() - run_start).totalNanoseconds() / 1e9
            assert expected_start - 1e-3 <= first_pulse < expected_start + 15.0
    finally:
        for name in set(mtd.getObjectNames()) - existing_workspaces:
            if mtd.doesExist(name):
                mtd.remove(name)


# CG2_9166: run end 300.048 s, last pulse at 300.048 s; the final slice ends one pulse period (1/60 s) later
_CG2_9166_EVENTS_END = 300.0649


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "time_interval, time_offset, expected_windows",
    [
        # list shorter than the run: the last value is repeated, the final slice absorbs the leftover
        ([100.0, 60.0], 10.0, [(10, 110), (110, 170), (170, 230), (230, _CG2_9166_EVENTS_END)]),
        # list longer than the run: the run end cuts the last slice short
        ([200.0, 150.0], 0.0, [(0, 200), (200, _CG2_9166_EVENTS_END)]),
        # single number: equal slices, the final slice absorbs the leftover (90.065 s)
        (70.0, 0.0, [(0, 70), (70, 140), (140, 210), (210, _CG2_9166_EVENTS_END)]),
    ],
    ids=["list shorter than the run", "list longer than the run", "single number"],
)
def test_time_slice_partition(datarepo_dir, time_interval, time_offset, expected_windows):
    """The run is partitioned with the last-slice rule, and each slice records its own window and length."""
    # Arrange: 300 s GPSANS run
    filename = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
    existing_workspaces = set(mtd.getObjectNames())

    try:
        # Act
        split_group = load_and_split(
            filename,
            sample_to_si_name="CG2:CS:SampleToSi",
            si_nominal_distance=0.0,
            time_interval=time_interval,
            time_offset=time_offset,
            output_workspace="time_slice_partition",
        )

        # Assert
        assert split_group.getNumberOfEntries() == len(expected_windows)
        for n, (start, end) in enumerate(expected_windows):
            sample_logs = SampleLogs(split_group.getItem(n))
            assert sample_logs.slice_start.value == pytest.approx(start, abs=1e-3)
            assert sample_logs.slice_end.value == pytest.approx(end, abs=1e-3)
            assert sample_logs.slice_interval.value == pytest.approx(end - start, abs=1e-3)
            assert sample_logs["slice_interval"].units == "seconds"
        if time_offset == 0.0:  # the slices hold every event of the run
            sliced_events = sum(split_group.getItem(n).getNumberEvents() for n in range(len(expected_windows)))
            assert sliced_events == 2313852
    finally:
        for name in set(mtd.getObjectNames()) - existing_workspaces:
            if mtd.doesExist(name):
                mtd.remove(name)


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "time_interval, time_offset, expected_rows",
    [(70.0, 0.0, 4), ([10.0, 20.0], 5.0, 15), ([0.7, 0.3], 0.0, 998), (0.1, 0.0, 3000), ([200.0, 150.0], 0.0, 2)],
    ids=["single number", "list with offset", "non-integer list", "short interval", "list longer than the run"],
)
def test_time_slice_splitter_rows(datarepo_dir, time_interval, time_offset, expected_rows):
    """The splitter table has one row per slice: Mantid never starts the list over for a last sliver of time."""
    # Arrange: 300 s GPSANS run
    filename = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
    existing_workspaces = set(mtd.getObjectNames())

    try:
        LoadEventNexus(Filename=filename, OutputWorkspace="time_slice_splitter_rows")
        time_filter = TimeIntervalFilter(
            "time_slice_splitter_rows", time_interval=time_interval, time_offset=time_offset
        )

        # Act
        time_filter.generate_filter()

        # Assert: the last row ends one pulse period after the last pulse
        splitter = mtd[time_filter.splitter_workspace]
        assert splitter.rowCount() == expected_rows
        run_start_ns = mtd["time_slice_splitter_rows"].run().startTime().totalNanoseconds()
        last_end = (splitter.cell(expected_rows - 1, 1) - run_start_ns) / 1e9
        assert last_end == pytest.approx(_CG2_9166_EVENTS_END, abs=1e-3)
    finally:
        for name in set(mtd.getObjectNames()) - existing_workspaces:
            if mtd.doesExist(name):
                mtd.remove(name)


@pytest.mark.datarepo
def test_time_slice_partition_with_load_time_cutoffs(datarepo_dir):
    """With a list and load-time cutoffs, the partition ignores the cutoffs and the empty slices are removed."""
    # Arrange: 300 s GPSANS run, slices of 20 s then 30 s; only events between 40 s and 100 s are loaded, so only
    # the slices 20-50 s, 50-80 s and 80-110 s hold events
    filename = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
    existing_workspaces = set(mtd.getObjectNames())

    try:
        # Act
        split_group = load_and_split(
            filename,
            sample_to_si_name="CG2:CS:SampleToSi",
            si_nominal_distance=0.0,
            time_interval=[20.0, 30.0],
            output_workspace="time_slice_partition_cutoffs",
            FilterByTimeStart=40.0,
            FilterByTimeStop=100.0,
        )
        cut_events = LoadEventNexus(
            Filename=filename,
            OutputWorkspace="time_slice_partition_cutoffs_reference",
            FilterByTimeStart=40.0,
            FilterByTimeStop=100.0,
        ).getNumberEvents()

        # Assert
        expected_windows = [(20, 50), (50, 80), (80, 110)]
        assert split_group.getNumberOfEntries() == len(expected_windows)
        for n, (start, end) in enumerate(expected_windows):
            sample_logs = SampleLogs(split_group.getItem(n))
            assert sample_logs.slice_start.value == pytest.approx(start)
            assert sample_logs.slice_end.value == pytest.approx(end)
            assert sample_logs.slice_interval.value == pytest.approx(30.0)
        # the kept slices hold every loaded event
        assert sum(split_group.getItem(n).getNumberEvents() for n in range(3)) == cut_events
    finally:
        for name in set(mtd.getObjectNames()) - existing_workspaces:
            if mtd.doesExist(name):
                mtd.remove(name)


if __name__ == "__main__":
    pytest.main([__file__])
