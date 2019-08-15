
from __future__ import print_function

import pytest

from mantid.simpleapi import (FindCenterOfMassPosition, LoadHFIRSANS,
                              MoveInstrumentComponent)
from ornl.sans.hfir.biosans import beam_finder


def test_beam_finder(biosans_f):
    '''
    Test with the new beam finder

    1. Find the beamcenter x,y
    2. Move detector1 x,y according to beamcenter x,y
    3. Find gravity
    4. Move wing_detector y according to Gravity

    '''

    ws = LoadHFIRSANS(Filename=biosans_f['beamcenter'])

    # 0.00144037741238 -0.0243732351545 -0.022044173994

    x, y, y_gravity = beam_finder.direct_beam_center(ws)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    print(x, y, y_gravity)
    assert x == pytest.approx(0.0014, abs=1e-3)
    assert y == pytest.approx(-0.0243, abs=1e-3)
    assert y_gravity == pytest.approx(-0.0220, abs=1e-3)

    # Let's center the instrument and get the new center:
    # Move down and left
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)

    # Now let's correct the wing detector for the gravity drop
    # Relative movement up words
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='wing_detector', X=0, Y=y_gravity)

    # After the re-centring we should be at (0,0)
    # Note that to give the same results we need to enter the center
    # estimates as the previous results!
    center = FindCenterOfMassPosition(InputWorkspace=ws,
                                      CenterX=-x, CenterY=-y)
    x, y = center
    # Tolerance 1e-3 == millimeters
    assert x == pytest.approx(0.0, abs=1e-3)
    assert y == pytest.approx(0.0, abs=1e-3)


def test_center_detector(biosans_f):

    ws = LoadHFIRSANS(Filename=biosans_f['beamcenter'])

    instrument = ws.getInstrument()
    pos_main = instrument.getComponentByName("detector1").getPos()
    pos_wing = instrument.getComponentByName("wing_detector").getPos()

    center_x = 0.0014
    center_y = -0.0243
    center_y_gravity = -0.0220

    ws = beam_finder.center_detector(ws, center_x, center_y, center_y_gravity)

    instrument = ws.getInstrument()
    pos_main_2 = instrument.getComponentByName("detector1").getPos()
    pos_wing_2 = instrument.getComponentByName("wing_detector").getPos()

    assert pytest.approx(abs(pos_main[0] - pos_main_2[0]), abs(center_x))
    assert pytest.approx(abs(pos_main[1] - pos_main_2[1]), abs(center_y))
    assert pytest.approx(
        abs(pos_wing[1] - pos_wing_2[1]), abs(center_y_gravity))
