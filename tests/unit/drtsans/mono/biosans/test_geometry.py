# local imports
from drtsans.mono.biosans.beam_finder import center_detector
from drtsans.geometry import get_position_south_detector
from drtsans.mono.biosans.geometry import (
    PHI_SPAN_MIDRANGE,
    set_position_south_detector,
    get_angle_wing_detector,
    set_angle_wing_detector,
    get_angle_midrange_detector,
    set_angle_midrange_detector,
    get_angle_south_detector,
    adjust_midrange_detector,
    midrange_to_wing_tubepixel,
)
from drtsans.mono.biosans.simulated_events import update_idf
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.simpleapi import LoadEmptyInstrument, LoadHFIRSANS, MoveInstrumentComponent
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# standard imports
#


def test_set_position_south_detector(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    # set the sample 0.042 meters away from the origin
    MoveInstrumentComponent(Workspace=workspace, ComponentName="sample-position", Z=-0.042, RelativePosition=False)
    set_position_south_detector(workspace, distance=7.0)  # meters
    assert_almost_equal(get_position_south_detector(workspace), 7.000, decimal=3)
    assert_almost_equal(SampleLogs(workspace).sample_detector_distance.value, 7.042, decimal=3)


def test_angle_wing_detector(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    assert_almost_equal(get_angle_wing_detector(workspace), 0.0, decimal=3)
    set_angle_wing_detector(workspace, angle=42.0)  # degrees
    assert_almost_equal(get_angle_wing_detector(workspace), 42.0, decimal=3)


def test_angle_midrange_detector(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    assert_almost_equal(get_angle_midrange_detector(workspace), 0.0, decimal=3)
    set_angle_midrange_detector(workspace, angle=42.0)  # degrees
    assert_almost_equal(get_angle_midrange_detector(workspace), 42.0, decimal=3)


def test_angle_south_detector(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    set_position_south_detector(workspace, distance=7.0)  # meters
    expected_angle = 4.3887  # arctan(0.537 / 7.0) and 0.537 is half the width of the south detector
    assert_almost_equal(get_angle_south_detector(workspace), expected_angle, decimal=3)


def test_midrange_rotation_by_criterium(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    set_position_south_detector(workspace, distance=7.0)  # meters
    # position the wing detector at an angle much bigger than the angle span of the midrange detector
    set_angle_wing_detector(workspace, angle=30.0)  # degrees
    phi = adjust_midrange_detector(workspace)
    assert_almost_equal(phi, get_angle_south_detector(workspace), decimal=3)
    # position the wing detector at an angle ensuring fair tube shadowing
    set_angle_wing_detector(workspace, get_angle_south_detector(workspace) + PHI_SPAN_MIDRANGE / 2)
    phi = adjust_midrange_detector(workspace)
    assert_almost_equal(phi, 4.0342, decimal=3)


def test_midrange_to_wing_tubepixel(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"))
    clean_workspace(workspace)
    SampleLogs(workspace).insert("wavelength", 12.0, unit="A")  # Angstroms
    md2ww = midrange_to_wing_tubepixel(workspace)
    assert [min(md2ww), max(md2ww)] == [90, 164]


def test_apply_samplelogs_midrange_rotation(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"))
    clean_workspace(workspace)
    # PV variable mr_rot_Readback follows the rotational convention of positive angles for clockwise rotation,
    # thus -42.0 means rotate the panel eastward (away from the beam) by 42 degrees
    SampleLogs(workspace).insert_time_series("mr_rot_Readback", [0.0], [-42.0], unit="deg")
    workspace = update_idf(workspace)
    assert_almost_equal(get_angle_midrange_detector(workspace), 42.0, decimal=2)


def test_apply_samplelogs_wing_rotation(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"))
    clean_workspace(workspace)
    SampleLogs(workspace).insert_time_series("ww_rot_Readback", [0.0], [42.0], unit="deg")
    workspace = update_idf(workspace)
    assert_almost_equal(get_angle_wing_detector(workspace), 42.0, decimal=2)


def test_apply_samplelogs_sample_detector_distance(fetch_idf, clean_workspace):
    workspace = LoadEmptyInstrument(Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"))
    clean_workspace(workspace)
    SampleLogs(workspace).insert_time_series("sample_detector_distance", [0.0], [12.0], unit="m")
    workspace = update_idf(workspace)
    assert_almost_equal(get_position_south_detector(workspace), 12.0, decimal=2)


@pytest.mark.datarepo
def test_api_geometry(biosans_f, clean_workspace):
    ws = LoadHFIRSANS(Filename=biosans_f["beamcenter"])
    clean_workspace(ws)

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
