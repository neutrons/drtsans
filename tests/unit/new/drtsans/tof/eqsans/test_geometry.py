import pytest
from pytest import approx
from drtsans.tof.eqsans.geometry import (sample_aperture_diameter,
                                           source_aperture_diameter,
                                           source_monitor_distance)
from drtsans.samplelogs import SampleLogs


def test_sample_aperture_diameter(serve_events_workspace):
    ws = serve_events_workspace('EQSANS_92353')
    sad = sample_aperture_diameter(ws)
    # ISSUE1887 TODO Enabled assert sad == approx(10)
    sad = SampleLogs(ws).single_value('sample-aperture-diameter')
    # ISSUE1887 TODO assert sad == approx(10)
    assert sad > 0


def test_source_aperture_diameter(serve_events_workspace):
    ws = serve_events_workspace('EQSANS_92353')
    sad = source_aperture_diameter(ws)
    # ISSUE187 TODO Enable assert sad == approx(20)
    sad = SampleLogs(ws).single_value('source-aperture-diameter')
    # ISSUE187 TODO Enable assert sad == approx(20)
    assert sad > 0


def test_source_monitor_distance(serve_events_workspace):
    ws = serve_events_workspace('EQSANS_92353')
    smd = source_monitor_distance(ws, unit='m')
    assert smd == approx(10.122, abs=0.001)
    smd = SampleLogs(ws).single_value('source-monitor-distance')
    assert smd == approx(10122, abs=1)


if __name__ == '__main__':
    pytest.main()
