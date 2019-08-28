import pytest
import os, numpy as np
from numpy.testing import assert_almost_equal
from mantid.simpleapi import LoadNexus
here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture(scope='module')
def testdata():
    path = os.path.join(here, 'Thickness_normalization_He.csv')
    data = np.genfromtxt(path, delimiter=',', skip_header=3)
    # intensity and errorbars, expected normalized intensity and errorbars
    I, E, normedI, normedE = data.T
    return I, E, normedI, normedE


@pytest.fixture(scope='module')
def histogram_workspace(testdata):
    I, E, normedI, normedE = testdata
    return


def test_thickness_normalization(histogram_workspace):
    return


if __name__ == '__main__':
    pytest.main()
