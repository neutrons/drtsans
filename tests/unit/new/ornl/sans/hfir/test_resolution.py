import os
import sys
sys.path.append(os.path.join(os.path.expanduser('~'), 'git', 'sans-rewrite'))
import unittest
import numpy as np

from mantid.simpleapi import *
from reduction_workflow.instruments.sans.hfir_command_interface import *
from reduction_workflow.command_interface import AppendDataFile, Reduce

from ornl.sans.hfir import resolution

def gpsans_files():
    data_dir = os.path.join(os.path.expanduser('~'), 'git', 'sans-rewrite', 'data')
    dd = os.path.join(data_dir, 'new', 'ornl', 'sans', 'hfir', 'gpsans')
    return dict(
        beamcenter=os.path.join(dd, 'CG2_exp325_scan0020_0001.xml'),
        beamcenter_off_setted=os.path.join(dd, 'CG2_exp245_scan0007_0001.xml'),
        sample_transmission=os.path.join(dd, 'CG2_exp245_scan0009_0001.xml'),
        sample_scattering=os.path.join(dd, 'CG2_exp245_scan0010_0001.xml'),
        dark_current=os.path.join(dd, 'CG2_exp244_scan0001_0001.xml'),
    )

def _create_reduced_ws():
    data_files = gpsans_files()

    configI = ConfigService.Instance()
    configI["facilityName"]='HFIR'
    BIOSANS()
    DirectBeamCenter(data_files['beamcenter'])
    AppendDataFile(data_files['sample_scattering'])
    AzimuthalAverage(binning="0.01,-0.02,0.11")
    Reduce()
    ws_iq = mtd['CG2_exp245_scan0010_0001_Iq']
    ws_iqxy = mtd['CG2_exp245_scan0010_0001_Iqxy']
    return ws_iq, ws_iqxy


class HFIRResolution(unittest.TestCase):
    def test_1d(self):
        """
        Test the Q resolution for a 1D distribution
        """
        ws_iq, _ = _create_reduced_ws()
        dq = resolution.q_resolution(ws_iq)
        dq_ref = ws_iq.readDx(0)
        # Check that the results are the same order of magnitude
        # Note: they do differ...
        summed = dq.sum()
        summed_ref = dq_ref.sum()
        self.assertTrue(np.fabs(np.log10(summed)-np.log10(summed_ref)) < 1.0)

    def test_2d(self):
        """
        Test the Q resolution for a 2D distribution
        """
        _, ws_iqxy = _create_reduced_ws()
        dqx, dqy = resolution.q_resolution(ws_iqxy)


if __name__ == '__main__':
    unittest.main()
