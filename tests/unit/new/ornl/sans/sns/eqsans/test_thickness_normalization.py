import pytest
import os, numpy as np
from numpy.testing import assert_almost_equal
from mantid.simpleapi import WorkspaceFactory
here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture(scope='module')
def testdata():
    path = os.path.join(here, 'Thickness_normalization_He.csv')
    data = np.genfromtxt(path, delimiter=',', skip_header=3)
    # intensity and errorbars, expected normalized intensity and errorbars
    I, E, normedI, normedE = data.T
    return I, E, normedI, normedE


@pytest.fixture(scope='module')
def test_workspaces(testdata):
    I, E, normedI, normedE = testdata
    nrows = 1
    nbins = I.size
    # input
    inputws = WorkspaceFactory.create(
        "Workspace2D", NVectors=nrows, XLength=nbins+1, YLength=nbins
    )
    inputws.setX(0, np.arange(nbins+1))
    inputws.setY(0, I)
    inputws.setE(0, E)
    # expected output
    expected_output_ws = WorkspaceFactory.create(
        "Workspace2D", NVectors=nrows, XLength=nbins+1, YLength=nbins
    )
    expected_output_ws.setX(0, np.arange(nbins+1))
    expected_output_ws.setY(0, normedI)
    expected_output_ws.setE(0, normedE)
    return inputws, expected_output_ws


def test_thickness_normalization(test_workspaces):
    inputws, expected_output_ws = test_workspaces
    print(expected_output_ws.dataY(0))
    return


if __name__ == '__main__':
    pytest.main()
