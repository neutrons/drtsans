# local imports
from drtsans.geometry import spectrum_info_ranges
from drtsans.instruments import empty_instrument_workspace
from drtsans.mono.biosans.geometry import (
    PIXELS_IN_TUBE,
    set_angle_wing_detector,
    set_angle_midrange_detector,
    set_position_south_detector,
)
from drtsans.samplelogs import SampleLogs
from drtsans.simulated_events import (
    insert_events,
    insert_beam_spot,
    insert_background,
    insert_events_isotropic,
    insert_events_ring,
    insert_events_sin_squared,
)


# third party imports
from mantid.api import mtd
from mantid.simpleapi import AddSampleLog, Integration
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# standard imports


def test_insert_events(temp_workspace_name):
    workspace = empty_instrument_workspace(
        temp_workspace_name(), filename="BIOSANS_Definition.xml", event_workspace=True
    )
    set_position_south_detector(workspace, distance=5.0)  # meters
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    AddSampleLog(
        Workspace=workspace,
        LogName="wavelength",
        LogText="18.0",
        LogType="Number Series",
        LogUnit="A",
        NumberType="Double",
    )

    def xy_intensities(x, y, _):
        return 1000

    insert_events(workspace, xy_intensities, components="detector1")
    first, next_to_last = spectrum_info_ranges(workspace, component="detector1")
    for index in range(first, next_to_last):
        assert workspace.getSpectrum(index).getNumberEvents() in [
            19,
            38,
            39,
            40,
        ]


def test_insert_beam_spot(temp_workspace_name):
    workspace = empty_instrument_workspace(
        temp_workspace_name(), filename="BIOSANS_Definition.xml", event_workspace=True
    )
    set_position_south_detector(workspace, distance=5.0)  # meters
    SampleLogs(workspace).insert("start_time", "2023-08-01 00:00:00")
    AddSampleLog(
        Workspace=workspace,
        LogName="wavelength",
        LogText="18.0",
        LogType="Number Series",
        LogUnit="A",
        NumberType="Double",
    )
    insert_beam_spot(workspace, center_x=-0.015, center_y=-0.03, diameter=0.01)
    assert_almost_equal(np.array(workspace.spectrumInfo().position(24952)), [-0.0141, -0.0306, 4.9958], decimal=3)
    assert workspace.getSpectrum(24952).getNumberEvents() == 39


def test_insert_background(temp_workspace_name):
    workspace_events = empty_instrument_workspace(
        temp_workspace_name(), filename="BIOSANS_Definition.xml", event_workspace=True
    )
    SampleLogs(workspace_events).insert("start_time", "2023-08-01 00:00:00")
    # Gaussian noise with mean 20 and stddev 4
    mean, stddev = 40, 4
    insert_background(workspace_events, flavor="gaussian noise", mean=mean, stddev=stddev)
    workspace_intensities = Integration(InputWorkspace=workspace_events, OutputWorkspace=temp_workspace_name())
    intensities = workspace_intensities.extractY().flatten()
    intensities_other = np.random.normal(mean, stddev, len(intensities)).astype(int)
    assert_almost_equal(np.mean(intensities), np.mean(intensities_other), decimal=1)
    assert_almost_equal(np.std(intensities), np.std(intensities_other), decimal=1)
    # Flat noise in between 2 and 7 counts
    workspace_events = empty_instrument_workspace(
        temp_workspace_name(), filename="BIOSANS_Definition.xml", event_workspace=True
    )
    SampleLogs(workspace_events).insert("start_time", "2023-08-01 00:00:00")
    min_counts, max_counts = 2, 7
    insert_background(workspace_events, flavor="flat noise", min_counts=min_counts, max_counts=max_counts)
    workspace_intensities = Integration(InputWorkspace=workspace_events, OutputWorkspace=temp_workspace_name())
    intensities = workspace_intensities.extractY().flatten()
    assert_almost_equal(np.min(intensities), min_counts, decimal=1)
    assert_almost_equal(np.max(intensities), max_counts, decimal=1)


@pytest.fixture(scope="function")
def biosans_workspace(temp_workspace_name, fetch_idf):
    workspace_events = empty_instrument_workspace(
        temp_workspace_name(), filename=fetch_idf("BIOSANS_Definition.xml"), event_workspace=True
    )
    SampleLogs(workspace_events).insert("start_time", "2023-08-01 00:00:00")
    AddSampleLog(
        Workspace=workspace_events,
        LogName="wavelength",
        LogText="18.0",
        LogType="Number Series",
        LogUnit="A",
        NumberType="Double",
    )
    set_position_south_detector(workspace_events, distance=5.0)  # meters
    return workspace_events


def _event_count_in_biosans_central_tube(input_workspace):
    r"""event count in tube 1 of bank13 (middle tube in the South detector)"""
    workspace = mtd[str(input_workspace)]
    first_id = 24576
    return sum([workspace.getSpectrum(first_id + i).getNumberEvents() for i in range(PIXELS_IN_TUBE)])


def test_insert_events_isotropic(biosans_workspace):
    insert_events_isotropic(
        biosans_workspace,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        component_efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_biosans_central_tube(biosans_workspace) == 868


def test_insert_events_ring(biosans_workspace):
    set_angle_midrange_detector(biosans_workspace, angle=5.0)  # degrees
    set_angle_wing_detector(biosans_workspace, angle=4.0)  # degrees
    insert_events_ring(
        biosans_workspace,
        twotheta_center=6.0,
        twotheta_dev=1.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        component_efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_biosans_central_tube(biosans_workspace) == 105


def test_insert_events_sin_squared(biosans_workspace):
    set_angle_midrange_detector(biosans_workspace, angle=5.0)  # degrees
    set_angle_wing_detector(biosans_workspace, angle=15.0)  # degrees
    insert_events_sin_squared(
        biosans_workspace,
        period=16.0,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        component_efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_biosans_central_tube(biosans_workspace) == 483


if __name__ == "__main__":
    pytest.main([__file__])
