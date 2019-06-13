#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadHFIRSANS
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.settings import unique_workspace_name

'''
import pytest
import os
os.chdir("/home/rhf/git/sans-rewrite")

pytest.main(["-vs", "/home/rhf/git/sans-rewrite/tests/unit/new/ornl/\
sans/hfir/gpsans/test_momentum_transfer.py"])
'''

def test_momentum_tranfer(gpsans_f):

    ws = LoadHFIRSANS(
        Filename=gpsans_f['sample_scattering'],
        OutputWorkspace=unique_workspace_name())

    wss_name_ws = bin_into_q2d(ws)
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)
