# local imports
from drtsans.instruments import (
    InstrumentEnumName,
    copy_to_newest_instrument,
    extract_run_number,
    empty_instrument_workspace,
    instrument_enum_name,
    is_time_of_flight,
)
from drtsans.mono.biosans.geometry import get_angle_wing_detector, get_position_south_detector
from drtsans.settings import amend_config, unique_workspace_dundername

# third party imports
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, LoadEmptyInstrument, LoadEventNexus
from numpy.testing import assert_almost_equal
import pytest

# standard imports
#


class TestInstrumentEnumName:
    def test_names(self):
        assert InstrumentEnumName.names() == ["BIOSANS", "EQSANS", "GPSANS"]


def test_instrument_name(serve_events_workspace):
    assert instrument_enum_name("EQ-SANS") == InstrumentEnumName.EQSANS
    assert str(instrument_enum_name("EQ-SANS")) == "EQSANS"
    assert str(instrument_enum_name("CG3")) == "BIOSANS"
    assert str(instrument_enum_name("somepath/CG3_961.nxs.h5")) == "BIOSANS"
    assert instrument_enum_name("nonexistantsansinstrument") == InstrumentEnumName.UNDEFINED
    assert instrument_enum_name(serve_events_workspace("EQSANS_92353.nxs.h5")) == InstrumentEnumName.EQSANS
    workspace = CreateWorkspace(DataX=range(42), DataY=range(42), OutputWorkspace=unique_workspace_dundername())
    assert instrument_enum_name(workspace) == InstrumentEnumName.UNDEFINED
    workspace.delete()


def test_is_time_of_flight(serve_events_workspace):
    for query in ("EQSANS", "EQ-SANS", serve_events_workspace("EQSANS_92353.nxs.h5")):
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
        output_workspace=unique_workspace_dundername(),
        instrument_name="BIOSANS",
        event_workspace=True,
        monitors_have_spectra=False,
    )
    detector_ids = workspace.detectorInfo().detectorIDs()
    detector_ids = detector_ids[detector_ids >= 0]
    assert workspace.getNumberHistograms() == len(detector_ids)
    assert list(workspace.getSpectrum(0).getDetectorIDs()) == [0]
    DeleteWorkspace(workspace)


def test_copy_to_newest_instrument(fetch_idf, reference_dir, clean_workspace):
    #
    # assert that spectra originally assigned to monitors is respected
    workspace1 = LoadEmptyInstrument(
        Filename=fetch_idf("BIOSANS_Definition_2019_2023.xml"), OutputWorkspace=unique_workspace_dundername()
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
        OutputWorkspace=unique_workspace_dundername(),
        MakeEventWorkspace=True,
    )
    clean_workspace(workspace2)
    old_histogram_count = workspace2.getNumberHistograms()
    workspace2 = copy_to_newest_instrument(workspace2)
    assert isinstance(workspace2, EventWorkspace)
    assert workspace2.getNumberHistograms() > old_histogram_count
    #
    # assert intensities and detector positions
    with amend_config(new_config={"default.instrument": "CG3"}, data_dir=reference_dir.new.biosans):
        workspace3 = LoadEventNexus(Filename="1322", OutputWorkspace=unique_workspace_dundername())
    clean_workspace(workspace3)
    pixel_counts = workspace3.getSpectrum(24956).getNumberEvents()
    position_detector1 = get_position_south_detector(workspace3)
    angle_wing = get_angle_wing_detector(workspace3)
    workspace3 = copy_to_newest_instrument(workspace3)
    assert workspace3.getSpectrum(24956).getNumberEvents() == pixel_counts
    assert_almost_equal(get_position_south_detector(workspace3), position_detector1, decimal=4)
    assert_almost_equal(get_angle_wing_detector(workspace3), angle_wing, decimal=3)


if __name__ == "__main__":
    pytest.main([__file__])
