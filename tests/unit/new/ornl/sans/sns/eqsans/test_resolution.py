"""
    Test EQSANS resolution. Requires Mantid nightly
"""
# import sys
# sys.path.insert(0, '/opt/mantidnightly/bin')  # noqa: E402
from __future__ import (absolute_import, division, print_function)

import os
import tempfile
import unittest

import numpy as np

from mantid import simpleapi as api
from ornl.sans.sns.eqsans import resolution
from reduction_workflow.command_interface import AppendDataFile, Reduce
from reduction_workflow.instruments.sans import sns_command_interface as eqsans


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
    eqsans.OutputPath(tempfile.gettempdir())
    Reduce()

    return api.mtd['EQSANS_92164.nxs']


class EQSANSResolution(unittest.TestCase):
    def test_2d(self):
        """
        Test the Q resolution for a 1D distribution
        """
        ws = _create_reduced_ws()
        dqx, dqy = resolution.q_resolution_per_pixel(ws)
        self.assertTrue(np.average(dqx) < 0.005)
        self.assertTrue(np.average(dqy) < 0.005)
        self.assertTrue(np.fabs(np.average(dqx) - np.average(dqy)) < 0.00002)


if __name__ == '__main__':
    unittest.main()
