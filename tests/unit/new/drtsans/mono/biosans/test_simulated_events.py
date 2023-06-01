# local imports
from drtsans.mono.biosans.simulated_events import (
    update_idf,
    insert_events_isotropic,
    insert_events_offset_center,
    insert_events_ring,
    insert_events_sin_squared,
)
from drtsans.mono.biosans.geometry import (
    set_angle_wing_detector,
    set_angle_midrange_detector,
    set_position_south_detector,
)
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.simpleapi import LoadEmptyInstrument, SaveNexus
import pytest

# standard imports


def test_update_idf(fetch_idf):
    idf_old = fetch_idf("BIOSANS_Definition_2019_2023.xml")
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=idf_old)
    assert workspace.getInstrument().getComponentByName("midrange_detector") is None
    workspace = update_idf(workspace)
    assert workspace.getInstrument().getComponentByName("midrange_detector")


def test_insert_events_isotropic(fetch_idf):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"), MakeEventWorkspace=True
    )
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector(workspace, distance=5.0)  # meters
    insert_events_isotropic(
        workspace,
        counts_per_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )


def test_insert_events_ring(fetch_idf):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"), MakeEventWorkspace=True
    )
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector(workspace, distance=5.0)  # meters
    set_angle_midrange_detector(workspace, angle=5.0)  # degrees
    set_angle_wing_detector(workspace, angle=4.0)  # degrees
    insert_events_ring(
        workspace,
        twotheta_center=6.0,
        twotheta_dev=1.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )


def test_insert_events_sin_squared(fetch_idf):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"), MakeEventWorkspace=True
    )
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector(workspace, distance=5.0)  # meters
    set_angle_midrange_detector(workspace, angle=5.0)  # degrees
    set_angle_wing_detector(workspace, angle=15.0)  # degrees
    insert_events_sin_squared(
        workspace,
        period=16.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        efficiencies=[1.0, 0.1, 1.0],
    )


def test_insert_events_offset_center(fetch_idf):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"), MakeEventWorkspace=True
    )
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    set_position_south_detector(workspace, distance=5.0)  # meters
    insert_events_offset_center(workspace, detector1_shift_x=-0.05, detector1_shift_y=-0.05)
    SaveNexus(InputWorkspace=workspace, Filename="/tmp/test_insert_events_offset_center.nxs")


if __name__ == "__main__":
    pytest.main([__file__])
