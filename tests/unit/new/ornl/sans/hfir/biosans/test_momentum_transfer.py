#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadHFIRSANS
from ornl.settings import unique_workspace_name
from ornl.sans.momentum_transfer import bin_into_q2d, bin_into_q1d


def test_momentum_tranfer(biosans_sensitivity_dataset):

    ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'],
        OutputWorkspace=unique_workspace_name())

    qxqy_wss_grouped = bin_into_q2d(ws, "detector1")
    assert len(qxqy_wss_grouped) == 3

    ws_iqxqy, ws_dqx, ws_dqy = qxqy_wss_grouped
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    qxqy_wss_grouped = bin_into_q2d(ws, "wing_detector")
    assert len(qxqy_wss_grouped) == 3

    ws_iqxqy, ws_dqx, ws_dqy = qxqy_wss_grouped
    assert ws_iqxqy.extractY().shape == (256, 160)
    assert ws_iqxqy.extractX().shape == (256, 161)

    ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)