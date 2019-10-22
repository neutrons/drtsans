import pytest
from pytest import approx
import drtsans.beam_finder as bf
from drtsans.settings import unique_workspace_dundername as uwd
from mantid.simpleapi import (CreateWorkspace, LoadInstrument, AddSampleLog,
                              MaskDetectors, LoadEmptyInstrument)
import numpy as np

# Note for testing beam center: The FindCenterOfMassPosition algorithm
# needs some pixels outside the given examples, to be able to perform
# the iterative approach. So, for a 2x2 example, we need a 4x4 instrument
# to padd with zeros on each side


@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4, 'dx': 1.,
                           'dy': 1., 'xc': 1.5, 'yc': 1.5}],
                         indirect=True)
def test_beam_finder_trivial(generic_IDF):
    r"""
    Testing section 3.1 in the master document
    Find beam center of a simple 2 x 2 detector
    Functions to test: drtsans.tof.eqsans.beam_finder.find_beam_center
    Underlying Mantid algorithms:
        FindCenterOfMassPosition https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
    dev - Andrei Savici <saviciat@ornl.gov>
    SME - Venky Pingali <pingalis@ornl.gov>
    """
    ws = CreateWorkspace(DataX=[1., 2.]*16,
                         DataY=[0., 0., 0., 0.,
                                0., 5., 3., 0.,
                                0., 11., 7., 0.,
                                0., 0., 0., 0.],
                         DataE=[1.]*16,
                         Nspec=16,
                         OutputWorkspace=uwd())
    LoadInstrument(Workspace=ws, InstrumentXML=generic_IDF,
                   RewriteSpectraMap=True, InstrumentName='EQ-SANS')
    AddSampleLog(Workspace=ws, LogName='wavelength', LogText='6.0', LogType='Number')
    inst = ws.getInstrument()
    assert inst.getDetector(5).getPos() == approx([2, 1, 5], abs=1e-5)
    assert inst.getDetector(6).getPos() == approx([2, 2, 5], abs=1e-5)
    assert inst.getDetector(9).getPos() == approx([1, 1, 5], abs=1e-5)
    assert inst.getDetector(10).getPos() == approx([1, 2, 5], abs=1e-5)
    x_cen, y_cen = bf.find_beam_center(ws)
    assert x_cen == approx(1.307692, abs=1e-5)
    assert y_cen == approx(1.384615, abs=1e-5)


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                           [0, 4, 2, 1, 3, 0, 2, 7, 1, 3, 1, 0],
                                           [0, 0, 4, 2, 3, 7, 5, 3, 4, 2, 2, 0],
                                           [0, 1, 2, 5, 9, 11, 9, 8, 6, 2, 3, 0],
                                           [0, 3, 4, 9, 20, 23, 28, 19, 9, 4, 1, 0],
                                           [0, 1, 3, 9, 27, 128, 79, 25, 13, 5, 1, 0],
                                           [0, 2, 4, 15, 23, 97, 201, 18, 15, 4, 2, 0],
                                           [0, 0, 1, 4, 18, 50, 41, 65, 5, 6, 2, 0],
                                           [0, 0, np.nan, np.nan, 7, 17, 16, 12, 8, 3, 3, 0],
                                           [0, 2,  np.nan, np.nan, 4, 3, 2, 4, 3, 2, 0, 0],
                                           [0, 0, 1, 1, 3, 2, 5, 1, 0, 1, 0, 0],
                                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                           'axis_values': [6., 7.],
                           'name':'BIOSANS', 'dx': 1, 'dy': 1, 'xc': 5.5, 'yc': 5.5}],
                         indirect=True)
def test_beam_finder_larger_workspace(generic_workspace):
    r"""
    Testing section 3.1 in the master document
    Find beam center of a larger 10(+2) x 10(+2) detector
    Functions to test: drtsans.tof.eqsans.beam_finder.find_beam_center
    Underlying Mantid algorithms:
        FindCenterOfMassPosition https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
    dev - Andrei Savici <saviciat@ornl.gov>
    SME - Venky Pingali <pingalis@ornl.gov>
    """
    ws = generic_workspace
    AddSampleLog(Workspace=ws, LogName='wavelength', LogText='6.0', LogType='Number')
    # masking
    intensities = ws.extractY()
    mask = np.where(np.isnan(intensities))[0].tolist()
    MaskDetectors(Workspace=ws, DetectorList=mask)
    # test masking
    inst = ws.getInstrument()
    assert inst.getDetector(mask[0]). getPos() == approx([3, 2, 5], abs=1e-5)
    assert inst.getDetector(mask[1]). getPos() == approx([3, 3, 5], abs=1e-5)
    assert inst.getDetector(mask[2]). getPos() == approx([2, 2, 5], abs=1e-5)
    assert inst.getDetector(mask[3]). getPos() == approx([2, 3, 5], abs=1e-5)
    spec = ws.spectrumInfo()
    for i in mask:
        assert spec.isMasked(i)
    # test functions
    x_cen, y_cen = bf.find_beam_center(ws)
    assert x_cen == approx(5.423913, abs=1e-5)
    assert y_cen == approx(5.654682, abs=1e-5)


def test_center_detector():
    r""" Testing moving EQSANS detector
    Functions to test: drtsans.tof.eqsans.beam_finder.center_detector
    Underlying Mantid algorithms:
        MoveInstrumentComponent https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html
    """
    # look at the original instrument
    w_eqsans = LoadEmptyInstrument(InstrumentName='EQ-SANS')
    inst = w_eqsans.getInstrument()
    assert inst.getDetector(0).getPos() == approx([0.524185, -0.520957, -0.0316256], abs=1e-5)
    assert inst.getDetector(49151).getPos() == approx([-0.525015, 0.520957, -0.0234656], abs=1e-5)

    # move detector
    xcenter = 0.5
    ycenter = 0.7
    bf.center_detector(w_eqsans, center_x=xcenter, center_y=ycenter)
    inst = w_eqsans.getInstrument()
    assert inst.getDetector(0).getPos() == approx([0.024185, -1.220957, -0.0316256], abs=1e-5)
    assert inst.getDetector(49151).getPos() == approx([-1.025015, -0.179043, -0.0234656], abs=1e-5)

# TODO test for moving wing detector


if __name__ == '__main__':
    pytest.main([__file__])
