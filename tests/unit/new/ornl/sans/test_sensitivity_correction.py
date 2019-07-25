import pytest
from mantid.kernel import Property
import numpy as np
from ornl.settings import unique_workspace_name as uwn
from ornl.sans.sensitivity import apply_sensitivity_correction,\
    calculate_sensitivity_correction
import os
from mantid.simpleapi import (mtd, DeleteWorkspace, LoadNexusProcessed)
import tempfile
from tests.conftest import data_dir

FILENAME = tempfile.NamedTemporaryFile('wt', suffix='.nxs').name
WKSPNAME_IN = 'EQSANS_87680'
FILENAME_IN = os.path.join(data_dir, 'EQSANS_87680_integrated.nxs')
MIN, MAX = 0.5, 2.0


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

    # check for the number masked
    info = wksp.detectorInfo()
    masked = [i for i in range(info.size())  # detectorInfo index != detId
              if info.isMasked(i)]
    assert len(masked) == 26513

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
        wksp = apply_sensitivity_correction(WKSPNAME_IN,
                                            sensitivity_filename=FILENAME,
                                            output_workspace=uwn())
    else:
        if not os.path.exists(FILENAME):
            raise RuntimeError('Sensitivity file "{}" '
                               'does not exist'.format(FILENAME))
        sensitivity = LoadNexusProcessed(Filename=FILENAME,
                                         OutputWorkspace=uwn())
        wksp = apply_sensitivity_correction(WKSPNAME_IN,
                                            sensitivity_workspace=sensitivity,
                                            output_workspace=uwn())
        sensitivity.delete()

    assert wksp

    origY = mtd[WKSPNAME_IN].extractY().flatten()
    newY = wksp.extractY().flatten()

    assert newY.size == origY.size
    assert not np.array_equal(newY, origY)


if __name__ == '__main__':
    pytest.main()
