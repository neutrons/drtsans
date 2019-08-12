
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

    x, y, y_gravity = beam_finder.direct_beam_center(ws)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
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
    center = FindCenterOfMassPosition(InputWorkspace=ws)
    x, y = center
    # Tolerance 1e-3 == milimeters
    assert x == pytest.approx(0.0, abs=1e-2)
    assert y == pytest.approx(0.0, abs=1e-3)
