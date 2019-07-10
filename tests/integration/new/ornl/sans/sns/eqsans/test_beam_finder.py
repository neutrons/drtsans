from tempfile import NamedTemporaryFile
import pytest
from pytest import approx
from mantid.simpleapi import (MoveInstrumentComponent, LoadEventNexus,
                              MaskBTP, ClearMaskFlag, ExtractMask, SaveMask)
from ornl.settings import amend_config, unique_workspace_dundername as uwd
from ornl.sans.sns import eqsans


def test_find_beam_center(eqsans_f, eqsans_p):
    r"""
    Integration test to find the location on the detector where
    the beam impinges

    1. Apply mask
    2. Find the beam center
    """
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        ws = LoadEventNexus(Filename=eqsans_f['beamcenter'],
                            OutputWorkspace=uwd())
    #
    # Find the beam center
    #
    assert eqsans.find_beam_center(ws) == approx((0.0265, 0.0180), abs=1e-04)
    #
    # Find the beam center with a mask workspace
    #
    eqsans.apply_mask(ws, Tube=eqsans_p['tubes_to_mask'])
    x0, y0 = eqsans.find_beam_center(ws)
    mask_ws = ExtractMask(ws, OutputWorkspace=uwd()).OutputWorkspace
    ClearMaskFlag(ws)
    assert eqsans.find_beam_center(ws, mask=mask_ws) == approx((x0, y0))
    #
    # Find the beam center with a mask file
    #
    with NamedTemporaryFile(delete=True, suffix='.xml') as f:
        SaveMask(InputWorkspace=mask_ws, OutputFile=f.name)
        xy = eqsans.find_beam_center(ws, mask=f.name)
        assert xy == approx((x0, y0), abs=1e-04)
    #
    # Let's move the beam center to the intersection point between the Z-axis
    # and the detector. The new (x, y) coordinates for the beam center
    # should be (0, 0) now.
    #
    eqsans.center_detector(ws, x=-x0, y=-y0)
    assert eqsans.find_beam_center(ws) == pytest.approx((0, 0), abs=1e-04)


if __name__ == '__main__':
    pytest.main()
