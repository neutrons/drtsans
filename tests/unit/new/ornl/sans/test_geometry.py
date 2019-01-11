from __future__ import (absolute_import, division, print_function)

import pytest
from mantid.simpleapi import LoadEmptyInstrument, MoveInstrumentComponent

from ornl.sans import geometry as geo


@pytest.fixture(scope='module')
def wss():
    r"""Just one workspace for each instrument"""

    # Load an EQSANS instrument and mess with the instrument components
    _eq_ws = LoadEmptyInstrument(InstrumentName='EQSANS')
    geo.sample_detector_distance(_eq_ws)
    for component, shift in (('detector1', 1.3),
                             ('sample-position', 0.02),
                             ('moderator', 0.3)):
        MoveInstrumentComponent(_eq_ws, ComponentName=component, Z=shift)

    # ssd: source-sample-distance, sdd: sample-detector-distance
    return dict(biosans=None,
                eqsans=dict(ws=_eq_ws, ssd=13842, sdd=1280),
                gpsans=None)


@pytest.mark.offline
def test_sample_source_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.sample_source_distance(v['ws']) ==\
                   pytest.approx(v['ssd'], rel=0.01)


@pytest.mark.offline
def test_sample_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.sample_detector_distance(v['ws']) ==\
                   pytest.approx(v['sdd'], rel=0.01)


@pytest.mark.offline
def test_source_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.source_detector_distance(v['ws']) ==\
                   pytest.approx(v['ssd'] + v['sdd'], rel=0.01)


if __name__ == '__main__':
    pytest.main()
