#!/usr/bin/env python
from __future__ import print_function

import tempfile
from os.path import join
from ornl.sans.save_ascii import save_ascii_1D, save_xml_1D

from mantid.simpleapi import Load
import numpy as np
import xml.etree.ElementTree as ET


def test_save_ascii_1d(reference_dir):
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel
    '''

    ws = Load(Filename=join(reference_dir.new.eqsans,
                            'test_save_output/EQSANS_68200_iq.nxs'))

    with tempfile.NamedTemporaryFile('wt') as tmp:
        save_ascii_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output = np.loadtxt(tmp.name, dtype={'names': ('Q', 'I', 'dI', 'dQ'),
                                             'formats': ('f', 'f', 'f', 'f')})
        reference = np.loadtxt(join(reference_dir.new.eqsans,
                                    'test_save_output/EQSANS_68200_Iq.txt'),
                               dtype={'names': ('Q', 'I', 'dI', 'dQ'),
                                      'formats': ('f', 'f', 'f', 'f')})
        assert np.allclose(output['Q'], reference['Q'], atol=1e-6)
        assert np.allclose(output['I'], reference['I'], atol=1e-6)
        assert np.allclose(output['dI'], reference['dI'], atol=1e-6)
        assert np.allclose(output['dQ'], reference['dQ'], atol=1e-6)

    with tempfile.NamedTemporaryFile('wt') as tmp:
        save_xml_1D(ws, 'test_reduction_log.hdf', tmp.name)
        output = []
        tree = ET.parse(tmp.name)
        root = tree.getroot()
        for node in root[0][2]:
            output.append([float(node[0].text), float(node[1].text),
                           float(node[2].text), float(node[3].text)])
        reference = []
        tree = ET.parse(join(reference_dir.new.eqsans,
                             'test_save_output/EQSANS_68200_Iq.xml'))
        root = tree.getroot()
        for node in root[0][2]:
            reference.append([float(node[0].text), float(node[1].text),
                              float(node[2].text), float(node[3].text)])
        assert np.allclose(output, reference, atol=1e-6)


def test_save_ascii_2d():
    pass
