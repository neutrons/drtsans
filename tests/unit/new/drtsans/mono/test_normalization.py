import pytest


@pytest.mark.offline
def test_normalisation_monitor(gpsans_f):

    from drtsans.mono.normalisation import monitor
    from drtsans.samplelogs import SampleLogs
    from mantid.simpleapi import LoadHFIRSANS
    from mantid import mtd

    input_sample_ws_mame = 'input_sample_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['sample_scattering'],
                 OutputWorkspace=input_sample_ws_mame)
    input_sample_ws = mtd[input_sample_ws_mame]
    # algorithm overwrites input workspace
    sample_logs = SampleLogs(input_sample_ws)
    monitor_counts = sample_logs.monitor.value
    test_value = input_sample_ws.readY(0)[0] * 10**8 / monitor_counts
    output_sample_ws = monitor(input_sample_ws)
    assert monitor_counts == 1284652
    assert output_sample_ws.readY(0)[0] == pytest.approx(test_value)


@pytest.mark.offline
def test_normalisation_time(gpsans_f):

    from drtsans.mono.normalisation import time
    from drtsans.samplelogs import SampleLogs
    from mantid.simpleapi import LoadHFIRSANS
    from mantid import mtd

    input_sample_ws_mame = 'input_sample_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['sample_scattering'],
                 OutputWorkspace=input_sample_ws_mame)
    input_sample_ws = mtd[input_sample_ws_mame]
    # algorithm overwrite input workspace
    sample_logs = SampleLogs(input_sample_ws)
    timer = float(sample_logs.timer.value)
    test_value = input_sample_ws.readY(612)[0] / timer
    output_sample_ws = time(input_sample_ws)
    assert timer == 60.0
    assert output_sample_ws.readY(612)[0] == pytest.approx(test_value)
