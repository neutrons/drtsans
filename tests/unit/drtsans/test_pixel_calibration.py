# local imports
from drtsans.pixel_calibration import loader_algorithm, BarPositionFormula, Table
from drtsans.settings import namedtuplefy

# third-party imports
from mantid.api import AnalysisDataService
import pytest

# standard imports
import json
import os
import tempfile
import shutil


@pytest.fixture(scope="session")
@namedtuplefy
def helper(reference_dir):
    r"""Helper object to access the database"""
    database_file = os.path.join(reference_dir.sans, "pixel_calibration", "calibrations.json")
    return {"database": database_file}


@pytest.fixture(scope="function")
def clone_database(helper, has_sns_mount):
    r"""Serve all contents of helper.database in a temporary database"""
    if not has_sns_mount:
        pytest.skip("Do not have /SNS properly mounted on this system")

    database_directory = os.path.dirname(helper.database)
    cloned_database_directory = tempfile.mkdtemp()
    os.rmdir(cloned_database_directory)  # shutil.copytree requires non-existing directory!
    shutil.copytree(database_directory, cloned_database_directory)
    cloned_database_file = os.path.join(cloned_database_directory, os.path.basename(helper.database))
    with open(cloned_database_file, "r") as file_handle:
        entries = json.load(file_handle)
        # replace the path of each table file
        for entry in entries:
            table_file = entry["tablefile"]
            assert os.path.exists(table_file)
            table_file = os.path.join(cloned_database_directory, "tables", os.path.basename(table_file))
            assert os.path.exists(table_file)
            entry["tablefile"] = table_file
    with open(cloned_database_file, "w") as file_handle:
        json.dump(entries, file_handle, indent=2)
    yield cloned_database_file
    # Tear down the temporary database
    shutil.rmtree(cloned_database_directory)


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "input_file, loader_name",
    [
        ("CG3_960.nxs.h5", "LoadEventNexus"),
        ("CG3_838.nxs", "LoadNexusProcessed"),
        ("BioSANS_exp327_scan0066_0001_mask.xml", "Load"),
    ],
)
def test_loader_algorithm(input_file, loader_name, datarepo_dir):
    input_file = os.path.join(
        datarepo_dir.biosans,
        "pixel_calibration",
        "test_loader_algorithm",
        input_file,
    )
    assert loader_algorithm(input_file).__name__ == loader_name


class TestBarPositionFormula:
    def test_elucidate_formula(self):
        formula = BarPositionFormula._elucidate_formula(("BIOSANS", "detector1"))
        assert formula == "565 - {y} + 0.0083115 * (191 - {tube})"
        with pytest.warns(UserWarning, match="Unable to find a bar position formula for argument Mary Poppings"):
            formula = BarPositionFormula._elucidate_formula("Mary Poppings")
        assert formula == "565 - {y} + 0.0 * {tube}"

    def test_validate_symbols(self):
        BarPositionFormula._validate_symbols("{y} {tube}")
        for invalid_formula in ("{tube}", "{dcal}", "y"):
            with pytest.raises(ValueError):
                BarPositionFormula._validate_symbols(invalid_formula)
        with pytest.warns(UserWarning, match='Formula does not contain "{tube}"'):
            assert BarPositionFormula._validate_symbols("565 - {y}") == "565 - {y} + 0.0 * {tube}"

    def test_str(self):
        with pytest.warns(UserWarning, match="Unable to find a bar position formula for argument unknown"):
            assert str(BarPositionFormula(instrument_component="unknown")) == BarPositionFormula._default_formula

    def test_evaluate(self):
        formula = BarPositionFormula(instrument_component=("GPSANS", "detector1"))  # use default formula
        assert formula.evaluate(565, 191) == pytest.approx(0.0)

    def test_validate_top_position(self):
        for formula in BarPositionFormula._default_formulae.values():
            BarPositionFormula(formula=formula).validate_top_position(0.0)
        with pytest.raises(RuntimeError):
            BarPositionFormula(("BIOSANS", "wing_detector")).validate_top_position(1150.0)
        with pytest.warns(UserWarning, match='Formula does not contain "{tube}"'):
            BarPositionFormula(formula="{y} - 565").validate_top_position(1150.0)


class TestTable:
    @pytest.mark.mount_eqsans
    def test_load(self, helper, has_sns_mount):
        r"""test method 'load'"""
        if not has_sns_mount:
            pytest.skip("Do not have /SNS properly mounted on this system")

        calibration = Table.load(helper.database, "BARSCAN", "GPSANS", "detector1", 20200104)
        assert calibration.daystamp == 20200103
        assert AnalysisDataService.doesExist("barscan_GPSANS_detector1_20200103")

    @pytest.mark.mount_eqsans
    def test_save(self, helper, clone_database, has_sns_mount):
        r"""test method 'save'"""
        if not has_sns_mount:
            pytest.skip("Do not have /SNS properly mounted on this system")

        calibration = Table.load(clone_database, "BARSCAN", "GPSANS", "detector1", 20200104)
        with pytest.raises(ValueError):
            calibration.save(database=clone_database)  # we cannot save a duplicate
        calibration.save(database=clone_database, overwrite=True)  # force saving a duplicate
        calibration = Table.load(clone_database, "BARSCAN", "GPSANS", "detector1", 20200104)
        assert os.path.dirname(calibration.tablefile) == os.path.join(os.path.dirname(clone_database), "tables")


if __name__ == "__main__":
    pytest.main([__file__])
