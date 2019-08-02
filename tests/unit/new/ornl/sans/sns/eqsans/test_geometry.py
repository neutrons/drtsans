import pytest
from pytest import approx
from mantid.simpleapi import LoadEventNexus, CloneWorkspace
from ornl.settings import amend_config, unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans.geometry import (sample_aperture_diameter,
                                           source_aperture_diameter,
                                           source_monitor_distance)
from ornl.sans.samplelogs import SampleLogs


@pytest.fixture(scope='function')
def ws(refd):
    if ws._original is None:
        with amend_config(data_dir=refd.new.eqsans):
            ws._original = LoadEventNexus('EQSANS_92353', OutputWorkspace=uwd())
            ws._all_ws.append(ws._original)
    clone = CloneWorkspace(ws._original, OutputWorkspace=uwd())
    ws._all_ws.append(clone)
    yield CloneWorkspace(ws._original, OutputWorkspace=uwd())
    [w.delete() for w in ws._all_ws]


ws._original, ws._all_ws = None, list()


def test_sample_aperture_diameter(ws):
    sad = sample_aperture_diameter(ws)
    assert sad == approx(10)
    sad = SampleLogs(ws).single_value('sample-aperture-diameter')
    assert sad == approx(10)


def test_source_aperture_diameter(ws):
    sad = source_aperture_diameter(ws)
    assert sad == approx(20)
    sad = SampleLogs(ws).single_value('source-aperture-diameter')
    assert sad == approx(20)


def test_source_monitor_distance(ws):
    smd = source_monitor_distance(ws, unit='m')
    assert smd == approx(10.122, abs=0.001)
    smd = SampleLogs(ws).single_value('source-monitor-distance')
    assert smd == approx(10122, abs=1)


if __name__ == '__main__':
    pytest.main()
