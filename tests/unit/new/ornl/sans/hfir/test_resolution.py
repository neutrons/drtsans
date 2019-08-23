"""
    Test EQSANS resolution. Requires Mantid nightly
"""
# import sys
# sys.path.insert(0, '/opt/mantidnightly/bin')  # noqa: E402
import os
import tempfile
import unittest

import numpy as np

from mantid import simpleapi as api
from ornl.sans.hfir import resolution
from ornl.sans.samplelogs import SampleLogs
from reduction_workflow.command_interface import AppendDataFile, Reduce
from reduction_workflow.instruments.sans import hfir_command_interface as hfir


def azimuthal_average(ws):
    """
        Test implementation of the azimuthal averaging so we can do
        a quick comparison of outputs.
    """
    sl = SampleLogs(ws)
    wl = sl.find_log_with_units('wavelength', 'Angstrom')

    spec_info = ws.spectrumInfo()

    q = np.zeros(spec_info.size())
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            q[i] = 4.0 * np.pi * np.sin(spec_info.twoTheta(i) / 2.0) / wl

    _, _, _dqx, _dqy = resolution.q_resolution_per_pixel(ws)
    bins = np.arange(0.01, 0.111, 0.001)
    iq, _ = np.histogram(q, bins=bins)
    dqx, _ = np.histogram(q, bins=bins, weights=_dqx)
    dqy, _ = np.histogram(q, bins=bins, weights=_dqy)

    dqx /= iq
    dqy /= iq
    dqx = np.nan_to_num(dqx)
    dqy = np.nan_to_num(dqy)
    dq = np.sqrt(dqx**2 + dqy**2)

    return dq


# TO-DO: use fixture reference_dir for accessing data files
def gpsans_files():
    from tests.conftest import data_dir
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
    hfir.AzimuthalAverage(binning="0.01,0.001,0.11")
    hfir.OutputPath(tempfile.gettempdir())
    Reduce()
    ws = api.mtd['CG2_exp245_scan0010_0001']
    ws_iq = api.mtd['CG2_exp245_scan0010_0001_Iq']
    ws_iqxy = api.mtd['CG2_exp245_scan0010_0001_Iqxy']
    return ws, ws_iq, ws_iqxy


class HFIRResolution(unittest.TestCase):

    def test_1d(self):
        """
        Test the Q resolution for a 1D distribution
        """
        _, ws_iq, _ = _create_reduced_ws()
        dq = resolution.q_resolution(ws_iq)
        dq_ref = ws_iq.readDx(0)
        self.assertTrue(dq.shape == dq_ref.shape == (100,))
        # Check that the results are the same order of magnitude
        # Note: they do differ...
        summed = dq.sum()
        summed_ref = dq_ref.sum()
        self.assertTrue(np.fabs(np.log10(summed)-np.log10(summed_ref)) < 1.0)

    def test_2d(self):
        """
        Test the Q resolution for a 2D distribution
        """
        _, _, ws_iqxy = _create_reduced_ws()
        dqx, dqy = resolution.q_resolution(ws_iqxy)
        self.assertTrue(dqx.shape == dqy.shape == (200, 200))
        self.assertTrue(np.average(dqx) < 0.03)
        self.assertTrue(np.average(dqy) < 0.03)

    def test_pixels(self):
        """
        Test the Q resolution for a 2D distribution
        """
        ws, _, _ = _create_reduced_ws()
        qx, qy, dqx, dqy = resolution.q_resolution_per_pixel(ws)

        self.assertTrue(qx.shape == qy.shape == dqx.shape == dqy.shape ==
                        (192*256+2,))
        self.assertTrue(np.min(np.abs(qx)) < np.max(np.abs(qx)))
        self.assertTrue(np.min(np.abs(qy)) < np.max(np.abs(qy)))

        self.assertTrue(np.average(dqx) < 0.03)
        self.assertTrue(np.average(dqy) < 0.03)

    def test_compare_Iq(self):
        """
            Test whether the averaged dq is similar to the reference.
        """
        # To ignore warning: "invalid value encountered in true_divide"
        np.seterr(divide='ignore', invalid='ignore')
        ws, ws_iq, _ = _create_reduced_ws()
        dq = resolution.q_resolution(ws_iq)
        dq_test = azimuthal_average(ws)
        self.assertTrue(np.fabs(np.average(dq) - np.average(dq_test)) < 0.002)


if __name__ == '__main__':
    unittest.main()
