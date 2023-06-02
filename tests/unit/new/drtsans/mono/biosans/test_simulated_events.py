# local imports
from drtsans.mono.biosans.simulated_events import (
    update_idf,
    insert_events_isotropic,
    insert_events_offset_center,
    insert_events_ring,
    insert_events_sin_squared,
)
from drtsans.mono.biosans.geometry import (
    PIXELS_IN_TUBE,
    set_angle_wing_detector,
    set_angle_midrange_detector,
    set_position_south_detector,
)
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.api import mtd
from mantid.simpleapi import DeleteWorkspace, LoadEmptyInstrument
import pytest

# standard imports


def _event_count_in_central_tube(input_workspace):
    r"""event count in tube 1 of bank13 (middle tube in the South detector)"""
    workspace = mtd[str(input_workspace)]
    first_id = 24276
    return sum([workspace.getSpectrum(first_id + i).getNumberEvents() for i in range(PIXELS_IN_TUBE)])


def test_update_idf(fetch_idf):
    idf_old = fetch_idf("BIOSANS_Definition_2019_2023.xml")
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=idf_old)
    assert workspace.getInstrument().getComponentByName("midrange_detector") is None
    workspace = update_idf(workspace)
    assert workspace.getInstrument().getComponentByName("midrange_detector")


def test_insert_events_isotropic(fetch_idf):
    LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace="test_insert_events_isotropic",
        MakeEventWorkspace=True,
    )
    SampleLogs("test_insert_events_isotropic").insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector("test_insert_events_isotropic", distance=5.0)  # meters
    insert_events_isotropic(
        "test_insert_events_isotropic",
        counts_per_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_central_tube("test_insert_events_isotropic") == 256
    DeleteWorkspace("test_insert_events_isotropic")


def test_insert_events_ring(fetch_idf):
    LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace="test_insert_events_ring",
        MakeEventWorkspace=True,
    )
    SampleLogs("test_insert_events_ring").insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector("test_insert_events_ring", distance=5.0)  # meters
    set_angle_midrange_detector("test_insert_events_ring", angle=5.0)  # degrees
    set_angle_wing_detector("test_insert_events_ring", angle=4.0)  # degrees
    insert_events_ring(
        "test_insert_events_ring",
        twotheta_center=6.0,
        twotheta_dev=1.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_central_tube("test_insert_events_ring") == 34
    DeleteWorkspace("test_insert_events_ring")


def test_insert_events_sin_squared(fetch_idf):
    LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace="test_insert_events_sin_squared",
        MakeEventWorkspace=True,
    )
    SampleLogs("test_insert_events_sin_squared").insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector("test_insert_events_sin_squared", distance=5.0)  # meters
    set_angle_midrange_detector("test_insert_events_sin_squared", angle=5.0)  # degrees
    set_angle_wing_detector("test_insert_events_sin_squared", angle=15.0)  # degrees
    insert_events_sin_squared(
        "test_insert_events_sin_squared",
        period=16.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_central_tube("test_insert_events_sin_squared") == 170
    DeleteWorkspace("test_insert_events_sin_squared")


def test_insert_events_offset_center(fetch_idf):
    LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace="test_insert_events_offset_center",
        MakeEventWorkspace=True,
    )
    SampleLogs("test_insert_events_offset_center").insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector("test_insert_events_offset_center", distance=5.0)  # meters
    insert_events_offset_center("test_insert_events_offset_center", detector1_shift_x=-0.05, detector1_shift_y=-0.05)
    assert _event_count_in_central_tube("test_insert_events_offset_center") == 4282
    DeleteWorkspace("test_insert_events_offset_center")


if __name__ == "__main__":
    pytest.main([__file__])
