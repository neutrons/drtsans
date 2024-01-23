# local imports
from drtsans.instruments import (
    InstrumentEnumName,
    copy_to_newest_instrument,
    extract_run_number,
    empty_instrument_workspace,
    instrument_enum_name,
    is_time_of_flight,
    instrument_facility_name,
)
from drtsans.instruments import fetch_idf as instruments_fetch_idf
from drtsans.mono.biosans.geometry import get_angle_wing_detector, get_position_south_detector

# third party imports
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, LoadEmptyInstrument, LoadEventNexus, mtd
from mantid.kernel import amend_config
from numpy.testing import assert_almost_equal
import pytest
from os.path import join as path_join

# standard imports
#


class TestInstrumentEnumName:
    def test_names(self):
        assert InstrumentEnumName.names() == ["BIOSANS", "EQSANS", "GPSANS"]


@pytest.mark.datarepo
def test_instrument_name(datarepo_dir, temp_workspace_name):
    assert instrument_enum_name("EQ-SANS") == InstrumentEnumName.EQSANS
    assert str(instrument_enum_name("EQ-SANS")) == "EQSANS"
    assert str(instrument_enum_name("CG3")) == "BIOSANS"
    assert str(instrument_enum_name("somepath/CG3_961.nxs.h5")) == "BIOSANS"
    assert instrument_enum_name("nonexistantsansinstrument") == InstrumentEnumName.UNDEFINED
    event_workspace = LoadEventNexus(
        path_join(datarepo_dir.eqsans, "EQSANS_92353.nxs.h5"), OutputWorkspace=temp_workspace_name()
    )
    assert instrument_enum_name(event_workspace) == InstrumentEnumName.EQSANS
    workspace = CreateWorkspace(DataX=range(42), DataY=range(42), OutputWorkspace=temp_workspace_name())
    assert instrument_enum_name(workspace) == InstrumentEnumName.UNDEFINED


@pytest.mark.datarepo
def test_is_time_of_flight(datarepo_dir, temp_workspace_name):
    event_workspace = LoadEventNexus(
        path_join(datarepo_dir.eqsans, "EQSANS_92353.nxs.h5"), OutputWorkspace=temp_workspace_name()
    )
    for query in ("EQSANS", "EQ-SANS", event_workspace):
        assert is_time_of_flight(query) is True
    for query in ("GPSANS", "CG2", "BIOSANS", "CG3"):
        assert is_time_of_flight(query) is False
    for query in [InstrumentEnumName.BIOSANS, InstrumentEnumName.GPSANS]:
        assert is_time_of_flight(query) is False
    for query in [InstrumentEnumName.EQSANS]:
        assert is_time_of_flight(query) is True


def test_extract_run_number():
    for query in (
        "/HFIR/CG3/CG3_961.nxs.h5",
        "CG3_961.nxs.h5",
        "CG3961",
        "CG3_961",
        "961",
        961,
    ):
        assert extract_run_number(query) == 961


def test_empty_instrument_workspace():
    workspace = empty_instrument_workspace(
        output_workspace=mtd.unique_hidden_name(),
        instrument_name="BIOSANS",
        event_workspace=True,
        monitors_have_spectra=False,
    )
    detector_ids = workspace.detectorInfo().detectorIDs()
    detector_ids = detector_ids[detector_ids >= 0]
    assert workspace.getNumberHistograms() == len(detector_ids)
    assert list(workspace.getSpectrum(0).getDetectorIDs()) == [0]
    DeleteWorkspace(workspace)


@pytest.mark.datarepo
def test_copy_to_newest_instrument(fetch_idf, datarepo_dir, clean_workspace):
    #
    # assert that spectra originally assigned to monitors is respected
    workspace1 = LoadEmptyInstrument(
        Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"), OutputWorkspace=mtd.unique_hidden_name()
    )
    clean_workspace(workspace1)
    old_histogram_count = workspace1.getNumberHistograms()
    workspace1 = copy_to_newest_instrument(workspace1)
    assert workspace1.getSpectrum(0).getDetectorIDs()[0] < 0
    assert workspace1.getNumberHistograms() > old_histogram_count
    #
    # assert that event workspace is preserved
    workspace2 = LoadEmptyInstrument(
        Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"),
        OutputWorkspace=mtd.unique_hidden_name(),
        MakeEventWorkspace=True,
    )
    clean_workspace(workspace2)
    old_histogram_count = workspace2.getNumberHistograms()
    workspace2 = copy_to_newest_instrument(workspace2)
    assert isinstance(workspace2, EventWorkspace)
    assert workspace2.getNumberHistograms() > old_histogram_count
    #
    # assert intensities and detector positions
    with amend_config(facility="HFIR", instrument="CG3", data_dir=datarepo_dir.biosans):
        workspace3 = LoadEventNexus(Filename="1322", OutputWorkspace=mtd.unique_hidden_name())
    clean_workspace(workspace3)
    pixel_counts = workspace3.getSpectrum(24956).getNumberEvents()
    position_detector1 = get_position_south_detector(workspace3)
    angle_wing = get_angle_wing_detector(workspace3)
    workspace3 = copy_to_newest_instrument(workspace3)
    assert workspace3.getSpectrum(24956).getNumberEvents() == pixel_counts
    assert_almost_equal(get_position_south_detector(workspace3), position_detector1, decimal=4)
    assert_almost_equal(get_angle_wing_detector(workspace3), angle_wing, decimal=3)


def test_fetch_idf(tmpdir):
    instruments_fetch_idf("BIOSANS_Definition_2019_2023.xml", output_directory=tmpdir)
    instruments_fetch_idf("BIOSANS_Definition.xml", output_directory=tmpdir)
    with pytest.raises(FileNotFoundError) as excinfo:
        instruments_fetch_idf("nonexisting.xml", output_directory=tmpdir)
    assert "nonexisting.xml" in str(excinfo.value)


@pytest.mark.datarepo
def test_instrument_facility_name(datarepo_dir, temp_workspace_name):
    # test instrument name
    assert instrument_facility_name("EQSANS") == "SNS"
    assert instrument_facility_name("BIOSANS") == "HFIR"
    assert instrument_facility_name("GPSANS") == "HFIR"
    assert instrument_facility_name("CG3") == "HFIR"
    assert instrument_facility_name("CG2") == "HFIR"
    with pytest.raises(ValueError):
        assert instrument_facility_name("")

    # test filename
    assert instrument_facility_name("somepath/CG3_961.nxs.h5") == "HFIR"
    assert instrument_facility_name("EQSANS_92353.nxs.h5") == "SNS"
    with pytest.raises(ValueError):
        assert instrument_facility_name("nonexistantsansinstrument")

    # test workspace
    event_workspace = LoadEventNexus(
        path_join(datarepo_dir.eqsans, "EQSANS_92353.nxs.h5"), OutputWorkspace=temp_workspace_name()
    )
    assert instrument_facility_name(event_workspace) == "SNS"
    workspace = CreateWorkspace(DataX=range(42), DataY=range(42), OutputWorkspace=temp_workspace_name())
    with pytest.raises(ValueError):
        assert instrument_facility_name(workspace)


if __name__ == "__main__":
    pytest.main([__file__])
