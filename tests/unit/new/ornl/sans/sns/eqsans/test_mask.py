import pytest
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadEmptyInstrument, ClearMaskFlag, MaskBTP,
                              CompareWorkspaces, ExtractMask,
                              SaveMask, LoadMask)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans.mask import apply_mask


def test_apply_mask():
    w = LoadEmptyInstrument(InstrumentName='EQ-SANS', OutputWorkspace=uwd())
    m = apply_mask(w, panel='front', Bank='25-48', Pixel='1-10',
                   output_workspace=uwd())
    assert isinstance(m, MaskWorkspace)
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, instrument='EQ-SANS', Components='front-panel')
    MaskBTP(Workspace=w, instrument='EQ-SANS', Bank='25-48', Pixel='1-10')
    m2 = ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result
    SaveMask(m, '/tmp/m.xml')
    SaveMask(m2, '/tmp/m2.xml')
    ClearMaskFlag(w)
    m = LoadMask(Instrument='EQ-SANS', InputFile='/tmp/m.xml',
                 RefWorkspace=w, OutputWorkspace='m')
    ClearMaskFlag(w)
    m2 = LoadMask(Instrument='EQ-SANS', InputFile='/tmp/m2.xml',
                  RefWorkspace=w, OutputWorkspace='m2')
    assert CompareWorkspaces(m, m2).Result


if __name__ == '__main__':
    pytest.main()
