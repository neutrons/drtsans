import os
import pytest
import tempfile
import shutil

from mantid.api import AnalysisDataService

from drtsans.pixel_calibration import Table
from drtsans.settings import namedtuplefy


@pytest.fixture(scope='session')
@namedtuplefy
def helper(reference_dir):
    database_file = os.path.join(reference_dir.new.sans, 'pixel_calibration', 'calibrations.json')
    return {'database': database_file}


@pytest.fixture(scope='function')
def clone_database(helper):
    r"""Serve all contents of helper.database in a temporary database"""
    database_directory = os.path.dirname(helper.database)
    cloned_database_directory = tempfile.mkdtemp()
    os.rmdir(cloned_database_directory)  # shutil.copytree requires non-existing directory!
    shutil.copytree(database_directory, cloned_database_directory)
    cloned_database_file = os.path.join(cloned_database_directory, os.path.basename(helper.database))
    yield cloned_database_file
    # Tear down the temporary database
    shutil.rmtree(cloned_database_directory)


class TestTable(object):

    def test_load(self, helper):
        r"""test method 'load'"""
        calibration = Table.load(helper.database, 'BARSCAN', 'GPSANS', 'detector1', 20200104)
        assert calibration.daystamp == 20200103
        assert AnalysisDataService.doesExist('barscan_GPSANS_detector1_20200103')

    def test_save(self, helper, clone_database):
        r"""test method 'save'"""
        calibration = Table.load(clone_database, 'BARSCAN', 'GPSANS', 'detector1', 20200104)
        assert os.path.dirname(calibration.tablefile) == os.path.join(os.path.dirname(helper.database), 'tables')
        with pytest.raises(ValueError):
            calibration.save(database=clone_database)  # we cannot save a duplicate
        calibration.save(database=clone_database, overwrite=True)  # force saving a duplicate
        calibration = Table.load(clone_database, 'BARSCAN', 'GPSANS', 'detector1', 20200104)
        assert os.path.dirname(calibration.tablefile) == os.path.join(os.path.dirname(clone_database), 'tables')


if __name__ == '__main__':
    pytest.main([__file__])