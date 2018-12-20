#!/usr/bin/env python
from __future__ import print_function

import pytest


def test_normalisation_monitor(gpsans_f):

    from ornl.sans.hfir.normalisation import monitor
    from ornl.sans.samplelogs import SampleLogs
    from mantid.simpleapi import LoadHFIRSANS
    from mantid import mtd

    input_sample_ws_mame = 'input_sample_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['sample_scattering'],
                OutputWorkspace=input_sample_ws_mame)
    input_sample_ws = mtd[input_sample_ws_mame]

    output_sample_ws = monitor(input_sample_ws)

    sample_logs = SampleLogs(input_sample_ws)
    monitor_counts = sample_logs.monitor.value
    assert monitor_counts == 1284652
    assert output_sample_ws.readY(0)[0] == \
        pytest.approx(input_sample_ws.readY(0)[0] / monitor_counts)


def test_normalisation_time(gpsans_f):

    from ornl.sans.hfir.normalisation import time
    from ornl.sans.samplelogs import SampleLogs
    from mantid.simpleapi import LoadHFIRSANS
    from mantid import mtd

    input_sample_ws_mame = 'input_sample_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['sample_scattering'],
                OutputWorkspace=input_sample_ws_mame)
    input_sample_ws = mtd[input_sample_ws_mame]

    output_sample_ws = time(input_sample_ws)

    sample_logs = SampleLogs(input_sample_ws)
    timer = float(sample_logs.timer.value)
    assert timer == 60.0
    assert output_sample_ws.readY(612)[0] == \
        pytest.approx(input_sample_ws.readY(612)[0] / timer)
