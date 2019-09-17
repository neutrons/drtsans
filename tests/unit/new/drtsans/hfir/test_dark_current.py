#!/usr/bin/env python
from __future__ import print_function
import pytest


@pytest.mark.offline
def test_dark_current(gpsans_f):
    from ornl.sans.hfir.normalisation import time
    from ornl.sans.hfir.dark_current import subtract_normalised_dark
    from ornl.sans.samplelogs import SampleLogs
    from mantid.simpleapi import LoadHFIRSANS
    from mantid import mtd

    # First read the data
    input_sample_ws_mame = 'input_sample_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['sample_scattering'],
                 OutputWorkspace=input_sample_ws_mame)
    sample_ws = mtd[input_sample_ws_mame]

    # second read and normalise the dark current
    input_dc_ws_mame = 'input_dc_ws_name'
    LoadHFIRSANS(Filename=gpsans_f['dark_current'],
                 OutputWorkspace=input_dc_ws_mame)
    input_dc_ws = mtd[input_dc_ws_mame]
    normalised_dc_ws = time(input_dc_ws)

    # third let's a DC subraction
    sample_subtracted = subtract_normalised_dark(
        sample_ws, normalised_dc_ws
    )

    # Let's test:
    normalised_sample_time = SampleLogs(sample_ws).timer.value

    v1_sample = sample_ws.dataY(612)[0]
    v1_dc = normalised_dc_ws.dataY(612)[0]
    v1_sample_subtracted = sample_subtracted.dataY(612)[0]

    assert v1_sample_subtracted == v1_sample - normalised_sample_time * v1_dc
