import pytest
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (
    LoadEmptyInstrument,
    ClearMaskFlag,
    MaskBTP,
    CompareWorkspaces,
    ExtractMask,
)
from drtsans.mask_utils import apply_mask


def test_apply_mask(clean_workspace, temp_workspace_name):
    w = LoadEmptyInstrument(InstrumentName="EQ-SANS", OutputWorkspace=temp_workspace_name())
    apply_mask(w, panel="front", Bank="25-48", Pixel="1-10")
    m = ExtractMask(w, OutputWorkspace=temp_workspace_name()).OutputWorkspace
    assert isinstance(m, MaskWorkspace)
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components="front-panel")
    MaskBTP(Workspace=w, Bank="25-48", Pixel="1-10")
    m2 = ExtractMask(w, OutputWorkspace=temp_workspace_name()).OutputWorkspace
    result_front, messages_front = CompareWorkspaces(m, m2)
    clean_workspace(messages_front)
    assert result_front
    #
    # Mask back panel
    #
    ClearMaskFlag(w)
    apply_mask(w, panel="back")
    m = ExtractMask(w, OutputWorkspace=temp_workspace_name()).OutputWorkspace
    ClearMaskFlag(w)
    MaskBTP(Workspace=w, Components="back-panel")
    m2 = ExtractMask(w, OutputWorkspace=temp_workspace_name()).OutputWorkspace
    result_back, messages_back = CompareWorkspaces(m, m2)
    clean_workspace(messages_back)
    assert result_back


if __name__ == "__main__":
    pytest.main([__file__])
