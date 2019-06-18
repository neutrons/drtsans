from __future__ import (absolute_import, division, print_function)

import pytest
from mantid.simpleapi import (MoveInstrumentComponent, LoadEventNexus,
                              MaskBTP)
from ornl.settings import amend_config, unique_workspace_name
from ornl.sans.sns.eqsans.beam_finder import direct_beam_center
from ornl.sans.sns import eqsans


def test_direct_beam_center(eqsans_f, eqsans_p):
    r"""
    Integration test to find the location on the detector where
    the beam impinges

    1. Apply mask
    2. Find the beam center
    """
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        ws = LoadEventNexus(Filename=eqsans_f['beamcenter'],
                            OutputWorkspace=unique_workspace_name())
    MaskBTP(Workspace=ws, Tube=eqsans_p['tubes_to_mask'])
    x, y = direct_beam_center(ws)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    assert (x, y) == pytest.approx((0.0265, 0.0180), abs=1e-04)
    # Let's center the instrument and get the new center: It should be
    # 0 after the re-centring
    MoveInstrumentComponent(ws, ComponentName='detector1', X=-x, Y=-y)
    assert direct_beam_center(ws) == pytest.approx((0, 0), abs=1e-04)


def test_direct_beam_from_top_level_api(eqsans_f, eqsans_p):
    """ Testing the top-level API """
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        ws = LoadEventNexus(Filename=eqsans_f['beamcenter'],
                            OutputWorkspace=unique_workspace_name())
    SANSMaskDTP(InputWorkspace=ws, tube=eqsans_p['tubes_to_mask'])
    x, y = eqsans.find_beam_center(ws)
    assert (x, y) == pytest.approx((0.0265, 0.0180), abs=1e-04)


if __name__ == '__main__':
    pytest.main()
