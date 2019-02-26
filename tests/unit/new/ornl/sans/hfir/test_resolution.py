import os
import unittest
import numpy as np
from mantid import simpleapi as api
from reduction_workflow.instruments.sans import hfir_command_interface as hfir
from reduction_workflow.command_interface import AppendDataFile, Reduce

from ornl.sans.hfir import resolution


def gpsans_files():
    _dir, _ = os.path.split(os.path.abspath(__file__))
    data_dir = os.path.join(_dir, '..', '..', '..', '..', '..', '..', 'data')
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

    configI = api.ConfigService.Instance()
    configI["facilityName"] = 'HFIR'
    hfir.BIOSANS()
    hfir.DirectBeamCenter(data_files['beamcenter'])
    AppendDataFile(data_files['sample_scattering'])
    hfir.AzimuthalAverage(binning="0.01,-0.02,0.11")
    Reduce()
    ws_iq = api.mtd['CG2_exp245_scan0010_0001_Iq']
    ws_iqxy = api.mtd['CG2_exp245_scan0010_0001_Iqxy']
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
        self.assertTrue(np.average(dqx) < 0.15)
        self.assertTrue(np.average(dqy) < 0.15)


if __name__ == '__main__':
    unittest.main()
