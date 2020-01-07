import numpy as np
import pytest

from mantid.simpleapi import LoadEmptyInstrument, MoveInstrumentComponent
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
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


@pytest.mark.parametrize('instrument, component, detmin, detmax',
                         [('EQ-SANS', '', 0, 49151),
                          ('BIOSANS', '', 0, 44 * 8 * 256 - 1),
                          ('BIOSANS', 'detector1', 0, 24 * 8 * 256 - 1),
                          ('BIOSANS', 'wing_detector', 24 * 8 * 256, 44 * 8 * 256 - 1)])
def test_bank_detector_ids(instrument, component, detmin, detmax):
    wksp = LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=unique_workspace_dundername())
    num_detectors = (detmax - detmin + 1)

    # None test
    detIDs = geo.bank_detector_ids(wksp, component=component, masked=None)
    assert detIDs.size == num_detectors
    assert detIDs.min() == detmin
    assert detIDs.max() == detmax

    detIDs = geo.bank_detector_ids(wksp, component=component, masked=False)
    assert detIDs.size == num_detectors
    assert detIDs.min() == detmin
    assert detIDs.max() == detmax

    detIDs = geo.bank_detector_ids(wksp, component=component, masked=True)
    assert len(detIDs) == 0


@pytest.mark.parametrize('instrument, component, wksp_index_min, wksp_index_max',
                         [('EQ-SANS', '', 1, 49151 + 2),
                          ('BIOSANS', '', 2, 44 * 8 * 256 + 2),
                          ('BIOSANS', 'detector1', 2, 24 * 8 * 256 + 2),
                          ('BIOSANS', 'wing_detector', 24 * 8 * 256 + 2, 44 * 8 * 256 + 2)])
def test_bank_workspace_indices(instrument, component, wksp_index_min, wksp_index_max):
    wksp = LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=unique_workspace_dundername())

    wksp_indices = geo.bank_workspace_index_range(wksp, component)
    assert wksp_indices[0] >= 0
    assert wksp_indices[1] <= wksp.getNumberHistograms()
    assert wksp_indices[0] == wksp_index_min
    assert wksp_indices[1] == wksp_index_max


def test_sample_aperture_diameter(serve_events_workspace):
    input_workspace = serve_events_workspace('EQSANS_92353')
    # diameter is retrieved from log 'beamslit4', and we convert the 10mm into 0.01 meters
    assert geo.sample_aperture_diameter(input_workspace) == pytest.approx(0.01, abs=0.1)
    # verify entry 'sample_aperture_diameter' has been added to the logs
    assert SampleLogs(input_workspace).single_value('sample_aperture_diameter') == pytest.approx(10.0, abs=0.1)


if __name__ == '__main__':
    pytest.main([__file__])
