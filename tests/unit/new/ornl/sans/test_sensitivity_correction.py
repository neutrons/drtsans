import pytest
from mantid.kernel import Property
import numpy as np
from ornl.settings import unique_workspace_name as uwn
from ornl.sans.sensitivity import apply_sensitivity_correction,\
    calculate_sensitivity_correction
import os
from mantid.simpleapi import (mtd, DeleteWorkspace, LoadNexusProcessed)
from tests.conftest import data_dir

FILENAME = os.path.join(os.path.abspath(os.curdir), 'test_sensitivity.nxs')
WKSPNAME_IN = 'EQSANS_87680'
FILENAME_IN = os.path.join(data_dir, 'EQSANS_87680_integrated.nxs')
MIN, MAX = 0.5, 2.0  # big because of bad choice of input file


@pytest.fixture(scope='module', autouse=True)
def eqsans_87680(request):
    '''Loads the data file and deletes it when all the tests are done'''
    LoadNexusProcessed(Filename=FILENAME_IN, OutputWorkspace=WKSPNAME_IN)

    yield  # return here

    # teardown code
    def fin():
        if WKSPNAME_IN in mtd:
            DeleteWorkspace(Workspace=WKSPNAME_IN)
    request.addfinalizer(fin)


def create_sensitivity_file():
    '''Single function to create the sensitivity file'''
    return calculate_sensitivity_correction(WKSPNAME_IN, MIN, MAX, FILENAME,
                                            output_workspace=uwn())


def test_calculation(cleanfile):
    wksp = create_sensitivity_file()
    cleanfile(FILENAME)
    assert wksp

    # verify resonable values in the output workspace
    values = wksp.extractY().flatten()
    values = values[values != Property.EMPTY_DBL]  # get rid of masked values
    assert values.size == 22639
    assert np.all(0.3 < values), 'Some values less than {}'.format(MIN)
    assert np.all(values < 1.5), 'Some values greater than {}'.format(MAX)

    # verify the output file exists
    assert os.path.exists(FILENAME)

    # cleanup
    wksp.delete()


@pytest.fixture(params=[True, False], ids=['from-file', 'from-workspace'])
def fromFile(request):
    return request.param


def test_apply_calculation(cleanfile, fromFile):
    if fromFile:
        if not os.path.exists(FILENAME):
            create_sensitivity_file()
            cleanfile(FILENAME)
        wksp = apply_sensitivity_correction(WKSPNAME_IN, filename=FILENAME,
                                            output_workspace=uwn())
    else:
        sensitivity = LoadNexusProcessed(Filename=FILENAME,
                                         OutputWorkspace=uwn())
        wksp = apply_sensitivity_correction(WKSPNAME_IN,
                                            sensitivity=sensitivity,
                                            output_workspace=uwn())
        sensitivity.delete()

    assert wksp

    origY = mtd[WKSPNAME_IN].extractY().flatten()
    newY = wksp.extractY().flatten()

    assert newY.size == origY.size
    assert not np.array_equal(newY, origY)


if __name__ == '__main__':
    pytest.main()
