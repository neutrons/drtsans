# local imports
from drtsans.mono.biosans.beam_finder import center_detector
from drtsans.mono.biosans.geometry import (
    PHI_SPAN_MIDRANGE,
    get_pixel_distances,
    get_position_south_detector,
    get_solid_angles,
    get_twothetas,
    set_position_south_detector,
    get_angle_wing_detector,
    set_angle_wing_detector,
    get_angle_midrange_detector,
    set_angle_midrange_detector,
    get_angle_south_detector,
    adjust_midrange_detector,
    midrange_to_wing_tubepixel,
)
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.simpleapi import LoadEmptyInstrument, LoadHFIRSANS, MoveInstrumentComponent
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# standard imports
#


def test_position_south_detector(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    # set the sample 0.042 meters away from the origing
    MoveInstrumentComponent(Workspace=workspace, ComponentName="sample-position", Z=-0.042, RelativePosition=False)
    set_position_south_detector(workspace, distance=7.0)  # meters
    assert_almost_equal(get_position_south_detector(workspace), 7.000, decimal=3)
    assert_almost_equal(SampleLogs(workspace).sample_detector_distance.value, 7.042, decimal=3)


def test_angle_wing_detector(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    assert_almost_equal(get_angle_wing_detector(workspace), 0.0, decimal=3)
    set_angle_wing_detector(workspace, angle=42.0)  # degrees
    assert_almost_equal(get_angle_wing_detector(workspace), 42.0, decimal=3)


def test_angle_midrange_detector(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    assert_almost_equal(get_angle_midrange_detector(workspace), 0.0, decimal=3)
    set_angle_midrange_detector(workspace, angle=42.0)  # degrees
    assert_almost_equal(get_angle_midrange_detector(workspace), 42.0, decimal=3)


def test_angle_south_detector(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    set_position_south_detector(workspace, distance=7.0)  # meters
    expected_angle = 4.3887  # arctan(0.537 / 7.0) and 0.537 is half the width of the south detector
    assert_almost_equal(get_angle_south_detector(workspace), expected_angle, decimal=3)


def test_midrange_rotation_by_criterium(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    set_position_south_detector(workspace, distance=7.0)  # meters
    # position the wing detector at an angle much bigger than the angle span of the midrange detector
    set_angle_wing_detector(workspace, angle=30.0)  # degrees
    phi = adjust_midrange_detector(workspace)
    assert_almost_equal(phi, get_angle_south_detector(workspace), decimal=3)
    # position the wing detector at an angle ensuring fair tube shadowing
    set_angle_wing_detector(workspace, get_angle_south_detector(workspace) + PHI_SPAN_MIDRANGE / 2)
    phi = adjust_midrange_detector(workspace)
    assert_almost_equal(phi, 4.0342, decimal=3)


def test_midrange_to_wing_tubepixel(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    SampleLogs(workspace).insert("wavelength", 12.0, unit="A")  # Angstroms
    md2ww = midrange_to_wing_tubepixel(workspace)
    assert [min(md2ww), max(md2ww)] == [90, 164]


def test_get_pixel_distances(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    set_position_south_detector(workspace, distance=7.0)  # meters
    expected_for_component = {
        "detector1": (7.0358, 7.0439),
        "wing_detector": (1.2402, 1.2477),
        "midrange_detector": (3.9962, 4.0043),
    }
    for component, expected in expected_for_component.items():
        distances = get_pixel_distances(workspace, component)
        assert_almost_equal((distances[0], distances[-1]), expected, decimal=3)


def test_get_solid_angles(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    set_position_south_detector(workspace, distance=7.0)  # meters
    expected_for_component = {
        "detector1": (0.02009, 0.0100),
        "wing_detector": (0.6097, 0.2994),
        "midrange_detector": (0.0626, 0.0311),
    }
    for component, expected in expected_for_component.items():
        solid_angles = get_solid_angles(workspace, component)
        assert_almost_equal((solid_angles[0], solid_angles[-1]), expected, decimal=3)
    #
    # move all components. Only solid angles for detector1 should change
    set_position_south_detector(workspace, distance=5.0)  # meters
    set_angle_wing_detector(workspace, angle=30.0)  # degrees
    set_angle_midrange_detector(workspace, angle=5.0)  # degrees
    expected_for_component["detector1"] = (0.03878, 0.0192)
    for component, expected in expected_for_component.items():
        solid_angles = get_solid_angles(workspace, component)
        assert_almost_equal((solid_angles[0], solid_angles[-1]), expected, decimal=3)


def test_get_twothetas(fetch_idf):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    set_position_south_detector(workspace, distance=7.0)  # meters
    expected_for_component = {
        "detector1": (6.105, 6.098),
        "wing_detector": (24.836, 49.498),
        "midrange_detector": (8.980, 7.475),
    }
    for component, expected in expected_for_component.items():
        twothetas = get_twothetas(workspace, component, units="degrees")
        assert_almost_equal((twothetas[0], twothetas[-1]), expected, decimal=2)


def test_api_geometry(biosans_f):
    ws = LoadHFIRSANS(Filename=biosans_f["beamcenter"])

    instrument = ws.getInstrument()
    pos_main = instrument.getComponentByName("detector1").getPos()
    pos_wing = instrument.getComponentByName("wing_detector").getPos()

    # Both detectors are at Y = 0
    assert pos_wing[1] == pos_main[1] == 0.0

    center_x = 0.0014
    center_y = -0.0243
    center_y_gravity = -0.0267

    center_detector(ws, center_x, center_y, center_y_gravity)

    instrument = ws.getInstrument()
    pos_main_2 = instrument.getComponentByName("detector1").getPos()
    pos_wing_2 = instrument.getComponentByName("wing_detector").getPos()

    np.testing.assert_allclose(abs(pos_main[0] - pos_main_2[0]), abs(center_x), rtol=1e-4)
    np.testing.assert_allclose(abs(pos_main[1] - pos_main_2[1]), abs(center_y), rtol=1e-4)
    np.testing.assert_allclose(abs(pos_wing[1] - pos_wing_2[1]), abs(center_y_gravity), rtol=1e-4)
    # Note that after the gravity correction the center Y of the wing detector
    # it's higher than the centre of the main detector
    assert pos_wing_2[1] > pos_main_2[1]


if __name__ == "__main__":
    pytest.main([__file__])
