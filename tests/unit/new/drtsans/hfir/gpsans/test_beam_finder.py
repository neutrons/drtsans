#!/usr/bin/env python
import pytest
from mantid import mtd
from mantid.simpleapi import MoveInstrumentComponent, FindCenterOfMassPosition, LoadHFIRSANS
from drtsans.settings import unique_workspace_dundername
from drtsans.hfir.gpsans import beam_finder


def test_beam_finder(gpsans_f):
    '''
    Test with the new beam finder
    '''

    ws_name = unique_workspace_dundername()
    LoadHFIRSANS(Filename=gpsans_f['beamcenter'], OutputWorkspace=ws_name)
    ws = mtd[ws_name]

    x, y = beam_finder.find_beam_center(ws)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))

    assert x == pytest.approx(0.02185, abs=1e-4)
    assert y == pytest.approx(-0.0193, abs=1e-4)

    # Let's center the instrument and get the new center: It should be 0 after
    # the re-centring
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)
    center = FindCenterOfMassPosition(InputWorkspace=ws)
    x, y = center

    # Tolerance 1e-3 == milimeters
    assert x == pytest.approx(0.0, abs=1e-4)
    assert y == pytest.approx(0.0, abs=1e-4)

    ws.delete()
