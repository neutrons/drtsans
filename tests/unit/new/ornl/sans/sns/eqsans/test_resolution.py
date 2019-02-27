import sys
import os
import unittest
import numpy as np
sys.path.insert(0, '/opt/mantidnightly/bin')  # noqa: E402
from mantid import simpleapi as api
from reduction_workflow.instruments.sans import sns_command_interface as eqsans
from reduction_workflow.command_interface import AppendDataFile, Reduce

from ornl.sans.sns.eqsans import resolution


def eqsans_files():
    ipts = '/SNS/EQSANS/IPTS-20196/nexus'
    # shared = '/SNS/EQSANS/shared/NeXusFiles/EQSANS'

    return dict(
        beamcenter=os.path.join(ipts, 'EQSANS_92160.nxs.h5'),
        sample_transmission=os.path.join(ipts, 'EQSANS_92162.nxs.h5'),
        sample_scattering=os.path.join(ipts, 'EQSANS_92164.nxs.h5'),
        dark_current=os.path.join(ipts, 'EQSANS_89157.nxs.h5'),
    )


def _create_reduced_ws():
    data_files = eqsans_files()

    configI = api.ConfigService.Instance()
    configI["facilityName"] = 'SNS'

    eqsans.EQSANS(False)
    AppendDataFile(data_files['sample_scattering'])
    eqsans.UseConfig(False)
    eqsans.UseConfigTOFTailsCutoff(False)
    eqsans.UseConfigMask(False)
    eqsans.SetBeamCenter(96.29, 126.15)
    eqsans.SetTransmission(1.0, 0.0)
    Reduce()

    return api.mtd['EQSANS_92164.nxs']


class EQSANSResolution(unittest.TestCase):
    def test_2d(self):
        """
        Test the Q resolution for a 1D distribution
        """
        ws = _create_reduced_ws()
        dqx, dqy = resolution.q_resolution_per_pixel(ws)
        self.assertTrue(np.average(dqx) < 0.006)
        self.assertTrue(np.fabs(np.average(dqx) -  np.average(dqy) < 0.00005))

if __name__ == '__main__':
    unittest.main()
