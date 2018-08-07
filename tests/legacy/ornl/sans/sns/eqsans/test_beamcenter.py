#!/usr/bin/env python

from __future__ import print_function

import sys
import os

import pytest

from dotenv import load_dotenv
load_dotenv()


def test_SANSBeamFinder():
    '''
    This is just a test for the legacy algorithm as it runs in mantid
    '''

    sys.path.append(os.getenv("MANTID_PATH"))
    from mantid.simpleapi import SANSBeamFinder
    from mantid.kernel import PropertyManagerDataService, PropertyManager

    pm = PropertyManager()
    PropertyManagerDataService.addOrReplace("test_pm", pm)

    out = SANSBeamFinder(
        Filename=os.path.join(os.getenv('DATA_DIRECTORY'),
                              'eqsans', 'EQSANS_68183_event.nxs'),
        # UseDirectBeamMethod=True,
        # BeamRadius=3,
        ReductionProperties='test_pm',
    )

    x = float(pm.getPropertyValue('LatestBeamCenterX'))
    y = float(pm.getPropertyValue('LatestBeamCenterY'))

    assert x == pytest.approx(90.6773456526)
    assert y == pytest.approx(131.698906123)

    assert out.FoundBeamCenterX == pytest.approx(90.6773456526)
    assert out.FoundBeamCenterY == pytest.approx(131.698906123)

    print(out.OutputMessage)


# todo
# FindCenterOfMassPosition
