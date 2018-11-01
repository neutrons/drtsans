#!/usr/bin/env python
from __future__ import print_function

import pytest


def test_beam_finder(eqsans_f, eqsans_p):
    '''
    Test with the new beam finder
    '''

    from ornl.sans.sns.eqsans import beam_finder
    from mantid import mtd
    from mantid.simpleapi import (
        MoveInstrumentComponent, FindCenterOfMassPosition)

    x, y = beam_finder.direct_beam_center(
        eqsans_f['beamcenter'], eqsans_p['tubes_to_mask'])
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    assert x == pytest.approx(0.02652545)
    assert y == pytest.approx(0.01804158)

    # Let's center the instrument and get the new center: It should be 0 after
    # the re-centring
    ws = mtd['ws_flattened']  # This comes from direct_beam_center function
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)
    center = FindCenterOfMassPosition(InputWorkspace=ws)
    x, y = center

    # Tolerance 1e-3 == milimeters
    assert x == pytest.approx(0.0, abs=1e-3)
    assert y == pytest.approx(0.0, abs=1e-3)
