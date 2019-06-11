from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pjn
from numpy.testing import assert_almost_equal

from mantid.kernel import ConfigService
from mantid.simpleapi import Load, EQSANSTofStructure, CompareWorkspaces


def test_eqsanstofstructure(refd):

    config = ConfigService.Instance()
    previous_instrument = config['instrumentName']
    config['instrumentName'] = 'EQSANS'
    previous_archive = config['datasearch.searcharchive']
    config['datasearch.searcharchive'] = 'hfir,sns'

    try:
        ddir = pjn(refd.legacy.eqsans, 'test_eqsanstofstructure')
        Load(Filename=pjn(ddir, 'EQSANS_92353_event_eqsanstofstr_in.nxs'),
             OutputWorkspace='change_tof_structure')
        alg_out = EQSANSTofStructure(InputWorkspace='change_tof_structure',
                                     LowTOFCut=500, HighTOFCut=2000)
        assert alg_out.FrameSkipping is True
        assert_almost_equal(alg_out.TofOffset, 11389.65, decimal=2)
        assert_almost_equal(alg_out.WavelengthMin, 2.6109, decimal=4)
        assert_almost_equal(alg_out.WavelengthMax, 5.7032, decimal=4)
        assert_almost_equal(alg_out.WavelengthMinFrame2, 9.5472, decimal=4)
        assert_almost_equal(alg_out.WavelengthMaxFrame2, 12.9809, decimal=4)
        Load(Filename=pjn(ddir, 'EQSANS_92353_event_eqsanstofstr_out.nxs'),
             OutputWorkspace='reference_workspace')
        cmp, mesg = CompareWorkspaces(Workspace1='change_tof_structure',
                                      Workspace2='reference_workspace',
                                      Tolerance=1e-4)
        assert cmp is True
    finally:
        config['instrumentName'] = previous_instrument
        config['datasearch.searcharchive'] = previous_archive


if __name__ == '__main__':
    pytest.main()
