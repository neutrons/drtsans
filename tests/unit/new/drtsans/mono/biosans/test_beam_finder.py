
from __future__ import print_function

import pytest

# https://docs.mantidproject.org/nightly/algorithms/LoadHFIRSANS-v1.html
# https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
# https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.htm
from mantid.simpleapi import (FindCenterOfMassPosition, LoadHFIRSANS,
                              MoveInstrumentComponent)
from drtsans.mono.biosans import beam_finder


def test_beam_finder(biosans_f):
    '''
    Test with the new beam finder

    1. Find the beamcenter x,y
    2. Move detector1 x,y according to beamcenter x,y
    3. Find gravity
    4. Move wing_detector y according to Gravity

    '''

    ws = LoadHFIRSANS(Filename=biosans_f['beamcenter'])

    # 0.00144037741238 -0.0243732351545 -0.0267
    x, y, y_gravity = beam_finder.find_beam_center(ws)

    assert x == pytest.approx(0.0014, abs=1e-3)
    assert y == pytest.approx(-0.0243, abs=1e-3)
    assert y_gravity == pytest.approx(-0.0267, abs=1e-3)

    instrument = ws.getInstrument()
    pos_main = instrument.getComponentByName("detector1").getPos()
    pos_wing = instrument.getComponentByName("wing_detector").getPos()

    # Let's center the instrument and get the new center:
    # Move down and left
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)

    pos_main_1 = instrument.getComponentByName("detector1").getPos()
    pos_wing_1 = instrument.getComponentByName("wing_detector").getPos()

    assert pos_main[0] - pos_main_1[0] == x
    assert pos_main[1] - pos_main_1[1] == y
    assert pos_wing == pos_wing_1  # we did not touch it

    # Now let's correct the wing detector for the gravity drop
    # Relative movement up words
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='wing_detector', X=-x, Y=-y_gravity)

    pos_main_2 = instrument.getComponentByName("detector1").getPos()
    pos_wing_2 = instrument.getComponentByName("wing_detector").getPos()

    assert pos_main_1 == pos_main_2
    assert pos_wing_2[1] == pos_main_2[1] + (abs(y_gravity) - abs(y))

    # After the re-centring we should be at (0,0)
    # Note that to give the same results we need to enter the center
    # estimates as the previous results!
    center = FindCenterOfMassPosition(InputWorkspace=ws,
                                      CenterX=-x, CenterY=-y)
    x1, y1 = center
    # Tolerance 1e-3 == millimeters
    assert x1 == pytest.approx(0.0, abs=1e-3)
    assert y1 == pytest.approx(0.0, abs=1e-3)

    # let's the test our wrap function. The results should be the same.
    x2, y2, y_gravity2 = beam_finder.find_beam_center(
        ws, CenterX=-x, CenterY=-y)

    assert x2 == pytest.approx(0.0, abs=1e-3) == x1
    assert y2 == pytest.approx(0.0, abs=1e-3) == y1
    assert y_gravity2 == pytest.approx(0.0 + y_gravity - y, abs=1e-3)
    assert abs(y_gravity2) > abs(y2)
