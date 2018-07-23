#!/usr/bin/env mantidpython

from __future__ import print_function

import pytest

from mantid.simpleapi import SANSBeamFinder
from mantid.kernel import PropertyManagerDataService, PropertyManager

pm = PropertyManager()
pmds = PropertyManagerDataService.add("test_pm", pm)


def test_default():
    '''
    This is just a test for the algorithm as it runs in mantid
    '''
    out = SANSBeamFinder(
        Filename='68183',
        # UseDirectBeamMethod=True,
        # BeamRadius=3,
        ReductionProperties='test_pm',
    )

    x = float(pm.getPropertyValue('LatestBeamCenterX'))
    y = float(pm.getPropertyValue('LatestBeamCenterY'))

    assert x == pytest.approx(90.6773456526)
    assert y == pytest.approx(131.698906123)

    assert out.FoundBeamCenterX == pytest.approx(90.6773456526)
    assert out.FoundBeamCenterX == pytest.approx(131.698906123)


# todo
# FindCenterOfMassPosition