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
def test_apply_mask_single_bin(generate_sans_generic_IDF):
    r"""Nine histograms, each containing only one bin"""
    wavelength = np.array([1, 2] * 9)
    intensities = np.array([8, 15, 30, 9, 17, 25, 6, 10, 31])
    ws = CreateWorkspace(DataX=wavelength,
                         DataY=intensities,
                         DataE=np.sqrt(intensities),
                         Nspec=9,
                         OutputWorkspace=uwd())
    LoadInstrument(Workspace=ws, InstrumentXML=generate_sans_generic_IDF,
                   RewriteSpectraMap=True, InstrumentName='GenericSANS')

    # detectors to be masked with detector ID's 0 and 3
    masked_detectors = [3, 0]
    apply_mask(ws, mask=masked_detectors)

    # check mask-flags
    si = ws.spectrumInfo()
    mask_flags = [False]*9
    for i in masked_detectors:
        mask_flags[i] = True
    assert [si.isMasked(i) for i in range(9)] == mask_flags

    # check masked data
    intensities_after_mask = np.copy(intensities)
    intensities_after_mask[masked_detectors] = 0
    assert ws.extractY().flatten() == approx(intensities_after_mask, abs=0.1)

    # errors also set to zero for the masked spectra
    errors_after_mask = np.sqrt(intensities_after_mask)
    assert ws.extractE().flatten() == approx(errors_after_mask, abs=1e-9)

    # check the value of a masked detector is irrelevant when doing
    # operations on the whole workspace
    big_value = 1.0e6
    ws.setY(3, [big_value])
    ws_sum = SumSpectra(ws, OutputWorkspace=uwd())
    assert ws_sum.dataY(0)[0] == sum(intensities_after_mask)
    intensities_after_mask[3] = big_value

    # Now apply an additional mask mimicking the trap beam
    x_c, y_c = (1, 1)  # location of the trap beam, in pixel coordinates
    pixels_per_tube = 3
    trap_detector_id = x_c * pixels_per_tube + y_c
    apply_mask(ws, mask=[trap_detector_id])
    intensities_after_mask[trap_detector_id] = 0
    assert ws.extractY().flatten() == approx(intensities_after_mask, abs=0.1)

    # Clean up
    ws.delete()
    ws_sum.delete()


@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 2, 'Ny': 2, 'dx': 0.00425,
                           'dy': 0.0055, 'xc': 0.32, 'yc': -0.16}],
                         indirect=True)
def test_apply_mask_simple_histogram(generate_sans_generic_IDF):
    r"""Four histograms, each containing three bins"""
    wavelength = np.array([1., 2., 3., 4.] * 4)
    intensities = np.array([[9, 10, 11, 3],
                            [8, 12, 4, 14],
                            [11, 15, 3, 16]]).transpose()
    ws = CreateWorkspace(DataX=wavelength, DataY=intensities,
                         DataE=np.sqrt(intensities),
                         Nspec=4, OutputWorkspace=uwd())
    LoadInstrument(Workspace=ws, InstrumentXML=generate_sans_generic_IDF,
                   RewriteSpectraMap=True, InstrumentName='GenericSANS')

    # detector to be masked with detector ID 2
    masked_detectors = [2]
    apply_mask(ws, mask=masked_detectors)

    # check mask-flags
    si = ws.spectrumInfo()
    mask_flags = [False, False, True, False]
    assert [si.isMasked(i) for i in range(4)] == mask_flags

    # check masked data
    intensities_after_mask = np.copy(intensities)
    intensities_after_mask[masked_detectors] = 0
    assert ws.extractY() == approx(intensities_after_mask, abs=0.1)

    # errors also set to zero for the masked spectra
    errors_after_mask = np.sqrt(intensities_after_mask)
    assert ws.extractE() == approx(errors_after_mask, abs=1e-6)

    # check the value of a masked detector is irrelevant when doing
    # operations on the whole workspace
    big_value = 1.0e6
    ws.setY(2, [big_value]*3)
    ws_sum = SumSpectra(ws, OutputWorkspace=uwd())
    assert ws_sum.dataY(0) == approx(np.sum(intensities_after_mask, axis=0))
    intensities_after_mask[2] = big_value

    # Now apply an additional mask mimicking the trap beam
    x_c, y_c = (1, 1)  # location of the trap beam, in pixel coordinates
    pixels_per_tube = 2
    trap_detector_id = x_c * pixels_per_tube + y_c
    apply_mask(ws, mask=[trap_detector_id])
    intensities_after_mask[trap_detector_id] = 0
    assert ws.extractY() == approx(intensities_after_mask, abs=0.1)

    # Clean up
    ws.delete()


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
