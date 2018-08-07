#!/usr/bin/env python
from __future__ import print_function

import sys
import os

import pytest

from dotenv import load_dotenv
load_dotenv()

sys.path.append(os.getenv("MANTID_PATH"))

# beam center file
FILENAME = os.path.abspath(
    os.path.join(
        os.getenv('DATA_DIRECTORY'), 'eqsans', 'EQSANS_68183_event.nxs'))
TUBES_TO_MASK = "1,48,53,54,85,123,130,137" #  From 1 to 192
    

def test_beam_finder():
    '''
    Test with the new beam finder
    '''

    from ornl.sans.sns.eqsans import beam_finder
    from mantid import mtd
    from mantid.simpleapi import (
        MoveInstrumentComponent, FindCenterOfMassPosition)

    x, y = beam_finder.direct_beam_center(FILENAME, TUBES_TO_MASK)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    assert x == pytest.approx(0.02652545)
    assert y == pytest.approx(0.01804158)

    # Let's center the instrument and get the new center: It should be 0 after
    # the re-centring
    ws = mtd['ws_flattened']  # This comes from direct_beam_center function
    MoveInstrumentComponent(Workspace=ws, ComponentName='detector1', X=-x, Y=-y)
    center = FindCenterOfMassPosition(InputWorkspace=ws)
    x, y = center

    # Tolerance 1e-3 == milimeters
    assert x == pytest.approx(0.0, abs=1e-3)
    assert y == pytest.approx(0.0, abs=1e-3)
