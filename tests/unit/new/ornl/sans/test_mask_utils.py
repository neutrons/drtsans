import pytest
from pytest import approx
import numpy as np
from tempfile import NamedTemporaryFile
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadEmptyInstrument, ClearMaskFlag, MaskBTP,
                              CompareWorkspaces, SaveMask, ExtractMask,
                              CreateWorkspace, LoadInstrument, SumSpectra)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.mask_utils import apply_mask


@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 3, 'Ny': 3, 'dx': 0.00425,
                           'dy': 0.0055, 'xc': 0.32, 'yc': -0.16}],
                         indirect=True)
def test_apply_mask_simple(generate_sans_generic_IDF):
    wavelength = np.array([1, 2] * 9)
    intensities = np.array([8, 15, 30, 9, 17, 25, 6, 10, 31])
    ws = CreateWorkspace(DataX=wavelength,
                         DataY=intensities,
                         DataE=np.sqrt(intensities),
                         Nspec=9)
    LoadInstrument(Workspace=ws, InstrumentXML=generate_sans_generic_IDF,
                   RewriteSpectraMap=True, InstrumentName='GenericSANS')
    # detectors to be masked with detector ID's 0 and 3
    masked_detectors = [3, 0]
    apply_mask(ws, mask=masked_detectors)
    data = ws.extractY()
    assert data == approx([0, 15, 30, 0, 17, 25, 6, 10, 31], abs=1e-9)
    err = np.sqrt(intensities)
    err[0] = 0
    err[3] = 0
    error = ws.extractE()
    assert error.flatten() == approx(err.flatten(), abs=1e-9)
    si = ws.spectrumInfo()
    for i in [3, 0]:
        assert si.isMasked(i) is True
    for i in range(9):
        if i not in [3, 0]:
            assert si.isMasked(i) is False
    ws.setY(3, [15000])  # the value of the masked detector is irrelevant
    ws_sum = SumSpectra(ws)
    assert ws_sum.dataY(0)[0] == sum(intensities) - intensities[0] - intensities[3]

    # apply an additional mask mimicking the trap beam
    x_c, y_c = (1, 1)  # location of the trap beam, in pixel coordinates
    trap_detector_id = x_c * 3 + y_c  # each tube has three detectors
    apply_mask(ws, mask=[trap_detector_id,])
    assert ws.extractY().flatten() == approx([0, 15, 30, 15000, 0, 25, 6, 10, 31], abs=1e-6)


@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 2, 'Ny': 2, 'dx': 0.00425,
                           'dy': 0.0055, 'xc': 0.32, 'yc': -0.16}],
                         indirect=True)
def test_apply_mask_simple(generate_sans_generic_IDF):
    wavelength = np.array([1., 2., 3., 4.] * 4)
    intensities = np.array([9, 10, 11, 3],
                           [8, 12, 4, 14],
                           [11, 15, 3, 16]).transpose()
    ws = CreateWorkspace(DataX=wavelength,
                         DataY=intensities,
                         DataE=np.sqrt(intensities),
                         Nspec=4)
    LoadInstrument(Workspace=ws, InstrumentXML=generate_sans_generic_IDF,
                   RewriteSpectraMap=True, InstrumentName='GenericSANS')
    masked_detectors = [2]
    apply_mask(ws, mask=masked_detectors)
    assert ws.readY(masked_detectors) == approx(np.zeros(3), abs=1e-6)

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
