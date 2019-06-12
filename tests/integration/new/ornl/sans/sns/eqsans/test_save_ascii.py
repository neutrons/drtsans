#!/usr/bin/env python
from __future__ import print_function

from mantid import mtd
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.settings import unique_workspace_name
from ornl.sans.save_ascii import save_ascii_1D, save_ascii_2D, save_xml_1D


'''
import pytest
import os
os.chdir("/home/rhf/git/sans-rewrite")

pytest.main(["-vs", "/home/rhf/git/sans-rewrite/tests/integration/new/ornl"\
    "/sans/sns/eqsans/test_momentum_transfer.py"])
'''


def bin_into_q2d_parallel(parameters):
    ws_name, component_name, out_ws_prefix = parameters
    ws = mtd[ws_name]  # need to pass the name. ws is shared between processes?
    workspaces = bin_into_q2d(ws, component_name, out_ws_prefix)
    return workspaces


def bin_into_q1d_parallel(parameters):
    ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix = parameters
    iq_ws = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix=out_ws_prefix)
    return iq_ws


# @pytest.mark.skip(reason="It only passes if run as standalone test")1
def test_save_ascii_parallel():
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel
    '''

    from mantid.simpleapi import AddSampleLog
    from ornl.sans.sns.eqsans import reduce
    ws = reduce.load_w('EQSANS_68200', unique_workspace_name(suffix="_raw"),
                       low_tof_clip=500, high_tof_clip=2000, dw=0.1)
    AddSampleLog(Workspace=ws, LogName='sample-aperture-diameter',
                 LogText='10.', LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=ws, LogName='source-aperture-diameter',
                 LogText='40.', LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=ws, LogName='wavelength',
                 LogText='4.8', LogType='Number', LogUnit='Angstrom')
    AddSampleLog(Workspace=ws, LogName='wavelength-spread',
                 LogText='0.7', LogType='Number', LogUnit='Angstrom')

    wss_name_ws = bin_into_q2d(ws)
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    save_ascii_1D(ws_iq, 'test_reduction_log.hdf',
                  '/tmp/test_1D_eqsans_lq.txt')
    with open('/tmp/test_1D_eqsans_lq.txt', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[101] == "0.246823\t0.000000\t0.000000\t0.035160\n"

    save_xml_1D(ws_iq, 'test_reduction_log.hdf', '/tmp/test_1D_eqsans_lq.xml')
    with open('/tmp/test_1D_eqsans_lq.xml', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[110] == '\t\t\t<Idata><Q unit="1/A">0.246823</Q><I '\
        'unit="Counts">0</I><Idev unit="Counts">0</Idev><Qdev unit='\
        '"1/A">0.0351595</Qdev></Idata>\n'

    save_ascii_2D(ws_iqxqy, ws_dqx, ws_dqy, 'test_reduction_log.hdf',
                  '/tmp/test_2D_eqsans_lq.dat')
    with open('/tmp/test_2D_eqsans_lq.dat', 'r') as output_file:
        output_lines = output_file.readlines()
    assert output_lines[48900] == "0.169809\t-0.177483\t0.000000\t0.000000\t"\
        "0.024772\t0.002118\n"
