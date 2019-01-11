#!/usr/bin/env python
from __future__ import print_function

import pytest


@pytest.mark.offline
def test_beam_finder(biosans_f):
    '''
    Test with the new beam finder
    '''

    from ornl.sans.hfir.biosans import beam_finder
    from mantid import mtd
    from mantid.simpleapi import (
        MoveInstrumentComponent, FindCenterOfMassPosition, LoadHFIRSANS)

    ws_name = "__beamcenter"
    LoadHFIRSANS(Filename=biosans_f['beamcenter'], OutputWorkspace=ws_name)
    ws = mtd[ws_name]

    x, y, y_gravity = beam_finder.direct_beam_center(ws)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    assert x == pytest.approx(0.0014, abs=1e-3)
    assert y == pytest.approx(-0.0243, abs=1e-3)
    assert y_gravity == pytest.approx(-0.0220, abs=1e-3)

    # Let's center the instrument and get the new center: It should be 0 after
    # the re-centring
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)
    center = FindCenterOfMassPosition(InputWorkspace=ws)
    x, y = center

    # Tolerance 1e-3 == milimeters
    assert x == pytest.approx(0.0, abs=1e-2)
    assert y == pytest.approx(0.0, abs=1e-3)
