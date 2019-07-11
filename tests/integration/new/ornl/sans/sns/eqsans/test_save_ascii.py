#!/usr/bin/env python
from __future__ import print_function

import tempfile
from os.path import join
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.sans.save_ascii import save_ascii_1D, save_ascii_2D, save_xml_1D

from mantid.simpleapi import Load
import numpy as np

def test_save_ascii_1d(refd):
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel
    '''

    ws = Load(Filename=join(refd.new.eqsans,'test_save_output/EQSANS_68200_iq.nxs'))

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output_lines = np.loadtxt(tmp.name,dtype={'names':('Q', 'I', 'dI', 'dQ'),'formats':('f', 'f', 'f', 'f')})
        reference_lines = np.loadtxt(join(refd.new.eqsans,'test_save_output/EQSANS_68200_Iq.txt'),
                                    dtype={'names':('Q', 'I', 'dI', 'dQ'),'formats':('f', 'f', 'f', 'f')})
        assert np.allclose(output_lines['Q'],reference_lines['Q'],atol=1e-6)
        assert np.allclose(output_lines['I'],reference_lines['I'],atol=1e-6)
        assert np.allclose(output_lines['dI'],reference_lines['dI'],atol=1e-6)
        assert np.allclose(output_lines['dQ'],reference_lines['dQ'],atol=1e-6)


    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_xml_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output_lines = tmp.readlines()
        assert output_lines[110] == '\t\t\t<Idata><Q unit="1/A">0.246823</Q>'\
            '<I unit="Counts">0</I><Idev unit="Counts">0</Idev><Qdev unit='\
            '"1/A">0.0351595</Qdev></Idata>\n'


def test_save_ascii_2d():
    pass
