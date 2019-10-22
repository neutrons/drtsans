import numpy as np
import pytest
from mantid.simpleapi import LoadEmptyInstrument, MoveInstrumentComponent
from drtsans.settings import unique_workspace_dundername
from drtsans import geometry as geo


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
def test_source_sample_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.source_sample_distance(v['ws']) ==\
                   pytest.approx(v['ssd'], rel=0.01)


@pytest.mark.offline
def test_sample_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.sample_detector_distance(v['ws']) == \
                   pytest.approx(v['sdd'], rel=0.01)


@pytest.mark.offline
def test_source_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.source_detector_distance(v['ws']) == \
                   pytest.approx(v['ssd'] + v['sdd'], rel=0.01)


def test_detector_translation():
    r"""Ascertain sub-components are moved when main detector is moved"""
    translation = np.array([0.01, 0.1, 1.0])
    detector_name = 'detector1'
    for instrument_name in ('EQ-SANS', 'CG2'):
        workspace = LoadEmptyInstrument(InstrumentName=instrument_name, OutputWorkspace=unique_workspace_dundername())
        instrument = workspace.getInstrument()
        component_detector = instrument.getComponentByName(detector_name)
        component_bank = instrument.getComponentByName('bank42')
        component_detector = instrument.getDetector(42)
        initial_positions = [c.getPos() for c in (component_detector, component_bank, component_detector)]
        MoveInstrumentComponent(workspace, ComponentName=detector_name,
                                RelativePosition=True, **dict(zip(('X', 'Y', 'Z'), translation)))
        final_positions = [c.getPos() for c in (component_detector, component_bank, component_detector)]
        for i, final_position in enumerate(final_positions):
            assert final_position == pytest.approx(np.array(initial_positions[i]) + translation, abs=1e-4)
        workspace.delete()


if __name__ == '__main__':
    pytest.main([__file__])
