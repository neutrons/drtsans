import pytest
from tempfile import NamedTemporaryFile
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadEmptyInstrument, ClearMaskFlag, MaskBTP,
                              CompareWorkspaces, SaveMask, ExtractMask)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.mask_utils import apply_mask


def test_apply_mask():
    w = LoadEmptyInstrument(InstrumentName='EQ-SANS', OutputWorkspace=uwd())
    #
    # using MaskBPT with a Component
    #
    m = apply_mask(w, Components='back-panel', output_workspace=uwd())
    assert isinstance(m, MaskWorkspace)
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components='back-panel')
    m2 = ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result
    ClearMaskFlag(w)
    m = apply_mask(w, Bank='1-24', Pixel='1-10', output_workspace=m.name())
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Bank='1-24', Pixel='1-10')
    m2 = ExtractMask(w, OutputWorkspace=m2.name()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result
    #
    # using a MaskWorkspace
    #
    ClearMaskFlag(w)
    m2 = apply_mask(w, mask=m, output_workspace=uwd())
    assert CompareWorkspaces(m, m2).Result
    m2 = ExtractMask(w, OutputWorkspace=m2.name()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result
    #
    # using a File
    #
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Bank='1-24', Pixel='1-10')
    m = ExtractMask(w, OutputWorkspace=m.name()).OutputWorkspace
    with NamedTemporaryFile(suffix='.xml') as f:
        SaveMask(m, f.name)
        m2 = apply_mask(w, mask=f.name, output_workspace=m2.name())
        assert CompareWorkspaces(m, m2).Result
    m2 = ExtractMask(w, OutputWorkspace=m2.name()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result
    #
    # using a MaskWorkspace and MaskBTP
    #
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components='front-panel')
    m2 = ExtractMask(w, OutputWorkspace=m2.name()).OutputWorkspace
    m = apply_mask(w, mask=m2, Bank='25-48', Pixel='1-10',
                   output_workspace=m.name())
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components='front-panel')
    MaskBTP(Workspace=w, Bank='25-48', Pixel='1-10')
    m2 = ExtractMask(w, OutputWorkspace=m2.name()).OutputWorkspace
    assert CompareWorkspaces(m, m2).Result


if __name__ == '__main__':
    pytest.main()