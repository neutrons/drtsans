import pytest
import os, numpy as np          # noqa: E401
from numpy.testing import assert_allclose
from mantid.simpleapi import WorkspaceFactory
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/thickness_normalization.py
from drtsans.thickness_normalization import normalize_by_thickness
here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture(scope='module')
def testdata():
    # create test data from csv file provided by LiLin He
    path = os.path.join(here, 'Thickness_normalization_He.csv')
    data = np.genfromtxt(path, delimiter=',', skip_header=3)
    # intensity and errorbars, expected normalized intensity and errorbars
    I, E, normedI, normedE = data.T
    return I, E, normedI, normedE


@pytest.fixture(scope='module')
def workspaces(testdata):
    """create workspace using the data loaded from the csv file provided by LiLin He"""
    # load data from csv
    I, E, normedI, normedE = testdata
    nrows = 1
    nbins = I.size
    # create input workspace
    inputws = WorkspaceFactory.create(
        "Workspace2D", NVectors=nrows, XLength=nbins+1, YLength=nbins
    )
    inputws.setX(0, np.arange(nbins+1))
    inputws.setY(0, I)
    inputws.setE(0, E)
    # create expected output workspace
    expected_output_ws = WorkspaceFactory.create(
        "Workspace2D", NVectors=nrows, XLength=nbins+1, YLength=nbins
    )
    expected_output_ws.setX(0, np.arange(nbins+1))
    expected_output_ws.setY(0, normedI)
    expected_output_ws.setE(0, normedE)
    return inputws, expected_output_ws


def test_thickness_normalization(workspaces):
    '''Test thickness normalization using data from a cvs file.
    The normalized result is compared to expected result provided in a csv file.

    Function tested: drtsans.tof.eqsans.api.normalize_by_thickness
    Underlying Mantid algorithms:
        Divide https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html

    dev - Jiao Lin <linjiao@ornl.gov>
    SME - LiLin He <hel3@ornl.gov>
    '''
    inputws, expected_output_ws = workspaces
    thickness = 0.1
    normed = normalize_by_thickness(inputws, thickness)
    assert_allclose(normed.readY(0), expected_output_ws.readY(0), rtol=5e-3)
    assert_allclose(normed.readE(0), expected_output_ws.readE(0), rtol=1e-7)
    return


if __name__ == '__main__':
    pytest.main([__file__])
