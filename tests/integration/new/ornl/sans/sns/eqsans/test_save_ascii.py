#!/usr/bin/env python
from __future__ import print_function

import tempfile
from os.path import join
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.sans.save_ascii import save_ascii_1D, save_ascii_2D, save_xml_1D

from mantid.simpleapi import Load

def test_save_ascii_1d(refd):
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel
    '''

    ws = Load(Filename=join(refd.new.eqsans,'test_save_output/eqsans_iq.nxs'))

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output_lines = tmp.readlines()
        assert output_lines[101] == "0.246823\t0.000000\t0.000000\t0.035160\n"

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_xml_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output_lines = tmp.readlines()
        assert output_lines[110] == '\t\t\t<Idata><Q unit="1/A">0.246823</Q>'\
            '<I unit="Counts">0</I><Idev unit="Counts">0</Idev><Qdev unit='\
            '"1/A">0.0351595</Qdev></Idata>\n'


def test_save_ascii_2d():
    pass
