# standard imports
from typing import Tuple, Union

import numpy.testing

# third party imports
from mantid.api import mtd
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import AddSampleLog, ConvertUnits, DeleteWorkspace, Integration, Rebin, SumSpectra
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# local imports
from drtsans.geometry import spectrum_info_ranges
from drtsans.instruments import empty_instrument_workspace, instrument_standard_name
from drtsans.mono.biosans.geometry import (
    PIXELS_IN_TUBE,
    set_angle_wing_detector,
    set_angle_midrange_detector,
    set_position_south_detector,
)
from drtsans.samplelogs import SampleLogs
from drtsans.settings import unique_workspace_dundername
from drtsans.simulated_events import (
    insert_events,
    insert_beam_spot,
    insert_background,
    insert_events_isotropic,
    insert_events_ring,
    insert_events_sin_squared,
)


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


@pytest.fixture(scope="function")
def eqsans_workspace(temp_workspace_name, fetch_idf):
    r"""
    Callable fixture returning an EventWorkspace with the latest EQSANS instrument definition

    Parameters
    ----------
    wavelength: float
        Log 'wavelength', to be used if the gravity drop correction is ticked on when generating events.
    start_time: str
        Log 'start_time' signaling the beginning of the synthetic run
    sample_detector_distance: float
        Position of the detector panels downstream the sample, in meters

    Returns
    -------
    EventWorkspace
    """

    def __eqsans_workspace(
        wavelength: float = 2.5, start_time: str = "2023-08-01 00:00:00", sample_detector_distance: float = 5.0
    ) -> EventWorkspace:
        workspace_events = empty_instrument_workspace(
            temp_workspace_name(), filename=fetch_idf("EQ-SANS_Definition.xml"), event_workspace=True
        )
        workspace_events.getAxis(0).setUnit("TOF")
        workspace_events.getAxis(1).setUnit("Label")
        SampleLogs(workspace_events).insert("start_time", start_time)
        AddSampleLog(
            Workspace=workspace_events,
            LogName="wavelength",
            LogText=f"{wavelength}",
            LogType="Number Series",
            LogUnit="A",
            NumberType="Double",
        )
        set_position_south_detector(workspace_events, distance=sample_detector_distance)  # meters
        return workspace_events

    return __eqsans_workspace


def _histogram_all_events(
    input_workspace: Union[str, EventWorkspace],
    units: str = "Wavelength",
    binning: str = "0.1,0.1,30",  # first bin boundary, bin width, last bin boundary
) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Convert events from time-of-flight to selected units, then histogram each spectrum and add all histograms, thus
    ending with just one histogram for all events.

    Parameters
    ----------
    input_workspace
    output_workspace
    units
    binning

    Returns
    -------
    Tuple with two items, representing the histogram's bin boundaries and the histogram intensities, respectively
    """
    temp_workspace = unique_workspace_dundername()  # a name not taken by any other already existing workspace
    ConvertUnits(InputWorkspace=input_workspace, OutputWorkspace=temp_workspace, Target=units)
    Rebin(InputWorkspace=temp_workspace, OutputWorkspace=temp_workspace, Params=binning, PreserveEvents=False)
    SumSpectra(InputWorkspace=temp_workspace, OutputWorkspace=temp_workspace)
    bins, intensities = mtd[temp_workspace].extractX()[0], mtd[temp_workspace].extractY()[0]
    DeleteWorkspace(temp_workspace)
    return bins, intensities


def _event_count_in_central_tube(input_workspace):
    r"""event count in tube 1 of bank13 (middle tube in the South detector)"""
    first_id = {"BIOSANS": 24576, "EQSANS": 24576}[instrument_standard_name(input_workspace)]
    return sum([mtd[str(input_workspace)].getSpectrum(first_id + i).getNumberEvents() for i in range(PIXELS_IN_TUBE)])


def test_insert_events_isotropic(biosans_workspace):
    insert_events_isotropic(
        biosans_workspace,
        max_counts_in_pixel=100,
        components=["detector1", "wing_detector", "midrange_detector"],
        component_efficiencies=[1.0, 0.1, 1.0],
    )
    assert _event_count_in_central_tube(biosans_workspace) == 868


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
    assert _event_count_in_central_tube(biosans_workspace) == 105


def test_insert_monocromatic_events_ring(eqsans_workspace):
    r"""
    Generate a set of neutron events imprinting a ring pattern on the detector. Verify that all events have the same
    neutron wavelength by converting from time-of-flight to wavelength and histogramming the result.
    """
    wavelength = 2.5  # units are Angstroms. All neutrons have this wavelength
    workspace = eqsans_workspace(wavelength=wavelength, sample_detector_distance=5.0)
    input_workspace = str(workspace)
    insert_events_ring(
        input_workspace,
        twotheta_center=3.0,
        twotheta_dev=0.5,
        max_counts_in_pixel=10,
        back_panel_attenuation=1.0,  # no attenuation
        solid_angle_correction=False,
        gravity_correction=False,
        components=["detector1"],
        component_efficiencies=[1.0],
        lambda_distribution=lambda events_count: np.repeat(wavelength, events_count),  # monochromatic distribution
    )
    numpy.testing.assert_allclose(_event_count_in_central_tube(input_workspace), 332, atol=10)

    # convert to wavelength and histogram between 0 and 5 Angstroms, with a bin width of 1.0 Angstrom. The result is
    # a histogram with five bins. The third bin has boundaries [2.0, 3.0] and therefore should contain all neutron
    # counts. The other bins should have zero counts.
    _, intensities = _histogram_all_events(input_workspace, units="Wavelength", binning="0.0,1.0,5.0")

    # indeterminacy below is atol + rtol * [0.0, 0.0, 0.0, 49728.0, 0.0]
    numpy.testing.assert_allclose(intensities, [0.0, 0.0, 49728.0, 0.0, 0.0], atol=0.1, rtol=0.01)


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
    assert _event_count_in_central_tube(biosans_workspace) == 483


if __name__ == "__main__":
    pytest.main([__file__])
