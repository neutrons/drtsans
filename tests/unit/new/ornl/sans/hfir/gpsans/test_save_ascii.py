#!/usr/bin/env python
from __future__ import print_function

import tempfile
from mantid.simpleapi import LoadHFIRSANS
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.settings import unique_workspace_name
from ornl.sans.save_ascii import save_ascii_1D, save_ascii_2D, save_xml_1D


def test_save_ascii(gpsans_f):

    ws = LoadHFIRSANS(
        Filename=gpsans_f['sample_transmission'],
        OutputWorkspace=unique_workspace_name())

    wss_name_ws = bin_into_q2d(ws)
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_1D(ws_iq, 'Test GPSANS', tmp.name)
        output_lines = tmp.readlines()
        print(output_lines[101])
        assert output_lines[101] == \
            "0.158594	0.000000	148.246656	0.025318\n"

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_xml_1D(ws_iq, 'Test GPSANS', tmp.name)
        output_lines = tmp.readlines()
        assert output_lines[110] == '\t\t\t<Idata><Q unit="1/A">0.158594</Q>'\
            '<I unit="Counts">0</I><Idev unit="Counts">0.235702</Idev>'\
            '<Qdev unit="1/A">0.0253178</Qdev></Idata>\n'

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_2D(ws_iqxqy, ws_dqx, ws_dqy, 'Test 2D GPSANS', tmp.name)
        output_lines = tmp.readlines()
        assert output_lines[48900] == "0.137103\t-0.081288\t0.000000\t"\
            "1.000000\t0.020047\t0.002688\n"
