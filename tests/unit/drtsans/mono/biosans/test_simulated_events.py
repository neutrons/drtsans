# local imports
from drtsans.mono.biosans.simulated_events import update_idf


# third party imports
from mantid.simpleapi import LoadEmptyInstrument
import pytest


def test_update_idf(fetch_idf, temp_workspace_name):
    idf_old = fetch_idf("BIOSANS_Definition_2019_2023.xml")
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=idf_old, OutputWorkspace=temp_workspace_name())
    assert workspace.getInstrument().getComponentByName("midrange_detector") is None
    workspace = update_idf(workspace)
    assert workspace.getInstrument().getComponentByName("midrange_detector")


if __name__ == "__main__":
    pytest.main([__file__])
