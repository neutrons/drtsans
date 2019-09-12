"""
    Test EQSANS resolution. Requires Mantid nightly
"""
# import sys
# sys.path.insert(0, '/opt/mantidnightly/bin')  # noqa: E402
from __future__ import (absolute_import, division, print_function)

import os
from os.path import join as pj
import tempfile
import unittest

import numpy as np

from mantid import simpleapi as api
from ornl.sans.sns.eqsans import momentum_transfer
from reduction_workflow.command_interface import AppendDataFile, Reduce
from reduction_workflow.instruments.sans import sns_command_interface as eqsans


def eqsans_files():
    from tests.conftest import data_dir
    dd = os.path.join(data_dir, 'new', 'ornl', 'sans', 'sns', 'eqsans')
    return dict(
        beamcenter=pj(dd, 'EQSANS_92160.nxs.h5'),
        sample_transmission=pj(dd, 'EQSANS_92162.nxs.h5'),
        sample_scattering=pj(dd, 'EQSANS_92164.nxs.h5'),
        dark_current=pj(dd, 'EQSANS_89157.nxs.h5'),
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
    eqsans.SetTransmission(1.0, 1.0)
    eqsans.OutputPath(tempfile.gettempdir())
    Reduce()

    return api.mtd['EQSANS_92164.nxs']


class EQSANSResolution(unittest.TestCase):
    def test_2d(self):
        """
        Test the Q and Q resolution for a 1D distribution
        """
        ws = _create_reduced_ws()
        qx, qy, dqx, dqy = momentum_transfer.q_resolution_per_pixel(ws)

        self.assertTrue(qx.shape == qy.shape == dqx.shape == dqy.shape)
        self.assertTrue(np.min(np.abs(qx)) < np.max(np.abs(qx)))
        self.assertTrue(np.min(np.abs(qy)) < np.max(np.abs(qy)))
        self.assertTrue(np.average(dqx) < 0.0055)
        self.assertTrue(np.average(dqy) < 0.0055)
        self.assertTrue(np.fabs(np.average(dqx) - np.average(dqy)) < 1e-4)

    def test_moderator_time_error(self):
        """Test moderator time uncertainty function using two wavelengths above and below 2 Angstroms
        and verify the output with expected results.
        dev - Jiao Lin <linjiao@ornl.gov>
        SME - William Heller <hellerwt@ornl.gov>

        For details see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/168
        """
        # the wavelengths to test
        wavelengths = [1.5, 9.3]
        # expected output
        expected = [214.74671875, 258.8954766]
        # calculate
        from ornl.sans.sns.eqsans.momentum_transfer import _moderator_time_error
        out = _moderator_time_error(np.array(wavelengths))
        # verify
        self.assertTrue(np.allclose(out, expected))
        return


if __name__ == '__main__':
    unittest.main()
