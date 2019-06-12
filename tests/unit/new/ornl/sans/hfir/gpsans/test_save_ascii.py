#!/usr/bin/env python
from __future__ import print_function

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

    save_ascii_1D(ws_iq, 'test_reduction_log.hdf',
                  '/tmp/test_1D_gpsans_lq.txt')
    with open('/tmp/test_1D_gpsans_lq.txt', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[101] == "0.158594\t0.000000\t0.235702\t0.025318\n"

    save_xml_1D(ws_iq, 'test_reduction_log.hdf', '/tmp/test_1D_gpsans_lq.xml')
    with open('/tmp/test_1D_gpsans_lq.xml', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[110] == '\t\t\t<Idata><Q unit="1/A">0.158594</Q><I '\
        'unit="Counts">0</I><Idev unit="Counts">0.235702</Idev><Qdev unit='\
        '"1/A">0.0253178</Qdev></Idata>\n'

    save_ascii_2D(ws_iqxqy, ws_dqx, ws_dqy, 'test_reduction_log.hdf',
                  '/tmp/test_2D_gpsans_lq.dat')
    with open('/tmp/test_2D_gpsans_lq.dat', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[48900] == "0.137103\t-0.081288\t0.000000\t1.000000\t"\
        "0.020047\t0.002688\n"
