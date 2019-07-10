import pytest
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadEmptyInstrument, ClearMaskFlag, MaskBTP,
                              CompareWorkspaces, ExtractMask)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans.mask import apply_mask


def test_apply_mask():
    w = LoadEmptyInstrument(InstrumentName='EQ-SANS', OutputWorkspace=uwd())
    m = apply_mask(w, panel='front', Bank='25-48', Pixel='1-10',
                   output_workspace=uwd())
    assert isinstance(m, MaskWorkspace)
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components='front-panel')
    MaskBTP(Workspace=w, Bank='25-48', Pixel='1-10')
    m2 = ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result


if __name__ == '__main__':
    pytest.main()
