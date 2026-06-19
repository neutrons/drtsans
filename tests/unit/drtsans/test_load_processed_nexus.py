"""
Unit tests for the allow_processed_nexus functionality in load_events.

This tests the ability to load processed Nexus files (saved with SaveNexusProcessed)
as a fallback when LoadEventNexus fails. This is used for live reduction where
the sample events are saved to a temporary file.
"""

import pytest
from mantid.simpleapi import (
    DeleteWorkspace,
    MoveInstrumentComponent,
    SaveNexusProcessed,
    mtd,
)
import numpy as np

from drtsans.instruments import empty_instrument_workspace
from drtsans.load import load_events
from drtsans.samplelogs import SampleLogs
from drtsans.simulated_events import insert_background


@pytest.fixture(scope="module")
def simulated_eqsans_events():
    """Create a simulated EQSANS instrument workspace with events.

    Returns
    -------
    EventWorkspace
        An EQSANS events workspace with simulated events
    """
    count = 9
    workspace_name = mtd.unique_hidden_name()
    workspace_events = empty_instrument_workspace(
        output_workspace=workspace_name, instrument_name="EQSANS", event_workspace=True
    )
    workspace_events.getAxis(0).setUnit("TOF")
    workspace_events.getAxis(1).setUnit("Label")
    MoveInstrumentComponent(Workspace=workspace_name, ComponentName="detector1", Z=5.0)
    sample_logs = SampleLogs(workspace_events)
    sample_logs.insert("start_time", "2023-08-01 00:00:00")
    sample_logs.insert("run_number", 99999)
    sample_logs.insert("experiment_identifier", "IPTS-99999")
    insert_background(
        workspace_events,
        # Normal distribution with 5 Angstroms mean wavelength, 0.1 Angstroms standard deviation
        lambda_distribution=lambda n_events: np.random.normal(loc=5.0, scale=0.1, size=n_events),
        flavor="fix count",
        flavor_kwargs={"count": count},
    )
    yield workspace_events
    # Cleanup
    if mtd.doesExist(workspace_name):
        DeleteWorkspace(workspace_name)


@pytest.fixture
def temp_processed_nexus_file(simulated_eqsans_events, tmp_path):
    """Save simulated events to a temporary processed Nexus file.

    Parameters
    ----------
    simulated_eqsans_events : EventWorkspace
        The simulated EQSANS events workspace
    tmp_path : Path
        Pytest fixture for temporary directory

    Returns
    -------
    str
        Path to the temporary processed Nexus file
    """
    temp_file = str(tmp_path / "EQSANS_99999_live.nxs")
    SaveNexusProcessed(InputWorkspace=simulated_eqsans_events, Filename=temp_file)
    return temp_file


def test_load_events_with_processed_nexus_fallback(temp_processed_nexus_file, temp_workspace_name):
    """Test that load_events can load a processed Nexus file when allow_processed_nexus=True.

    This simulates the live reduction scenario where events are saved to a temporary
    processed Nexus file.
    """
    ws_name = temp_workspace_name()

    # Load the processed Nexus file using the fallback mechanism
    ws = load_events(
        run=temp_processed_nexus_file,
        output_workspace=ws_name,
        allow_processed_nexus=True,
    )

    # Verify the workspace was loaded successfully
    assert ws is not None
    assert mtd.doesExist(ws_name)
    # Check basic workspace properties
    assert ws.getNumberHistograms() > 0


def test_load_events_without_processed_nexus_fails(temp_processed_nexus_file, temp_workspace_name):
    """Test that load_events fails when trying to load a processed Nexus file
    without allow_processed_nexus=True.

    This is the expected behavior for normal autoreduction which uses event Nexus files.
    """
    ws_name = temp_workspace_name()

    # Attempting to load without allow_processed_nexus should fail
    with pytest.raises(RuntimeError):
        load_events(
            run=temp_processed_nexus_file,
            output_workspace=ws_name,
            allow_processed_nexus=False,
        )


def test_load_events_allow_processed_nexus_preserves_workspace_properties(
    simulated_eqsans_events, temp_processed_nexus_file, temp_workspace_name
):
    """Test that loading a processed Nexus file preserves workspace properties.

    Verifies that SaveNexusProcessed + Load gives back an equivalent workspace.
    """
    ws_name = temp_workspace_name()

    # Load the processed Nexus file
    loaded_ws = load_events(
        run=temp_processed_nexus_file,
        output_workspace=ws_name,
        allow_processed_nexus=True,
    )

    # Compare with the original workspace
    original_ws = simulated_eqsans_events

    # Check that the number of histograms matches
    assert loaded_ws.getNumberHistograms() == original_ws.getNumberHistograms()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
