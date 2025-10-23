import os
from unittest.mock import Mock

from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import DeleteWorkspace, MoveInstrumentComponent, mtd, Rebin
import numpy as np
import pytest

from .script_locator import reduce_EQSANS
from drtsans.instruments import empty_instrument_workspace
from drtsans.samplelogs import SampleLogs
from drtsans.simulated_events import insert_background


@pytest.fixture(scope="module")
def simulated_events() -> EventWorkspace:
    """Create a simulated EQSANS instrument workspace with 9 events per pixel

    Also adds sample logs:
    - start_time = 2023-08-01 00:00:00
    - run_number = 12345

    Parameters
    ----------
    temp_workspace_name : callable
        A fixture that returns a unique workspace name when called

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
    sample_logs.insert("run_number", 12345)
    sample_logs.insert("experiment_identifier", "IPTS-12345")
    insert_background(
        workspace_events,
        # Normal distribution with 5 Angstroms mean wavelength, 0.1 Angstroms standard deviation
        lambda_distribution=lambda n_events: np.random.normal(loc=5.0, scale=0.1, size=n_events),
        flavor="fix count",
        flavor_kwargs={"count": count},
    )
    workspace_events = Rebin(
        InputWorkspace=workspace_name,
        Params=[0.0, 1000.0, 50000.0],
        OutputWorkspace=workspace_name,
    )
    yield workspace_events
    # Cleanup
    DeleteWorkspace(workspace_name)


@pytest.mark.mount_eqsans
def test_autoreduce_sample(simulated_events, tmp_path):
    mock_args = Mock()
    mock_args.events_file = "/SNS/EQSANS/IPTS-34577/nexus/EQSANS_162568.nxs.h5"
    mock_args.outdir = str(tmp_path)
    mock_args.no_publish = True  # Disable publishing to the live data server
    reduce_EQSANS.autoreduce(mock_args)
    filenames = [
        "autoreduce_162568.log",
        "EQSANS_162568.html",
        "EQSANS_162568_Iq.dat",
        "EQSANS_162568_Iq.png",
        "EQSANS_162568_Iqxqy.dat",
        "EQSANS_162568_Iqxqy.h5",
        "EQSANS_162568_Iqxqy.png",
        "EQSANS_162568_processed.nxs",
        "EQSANS_162568_reduction_log.hdf",
        "reduction_options_162568.json",
    ]
    for expected in filenames:
        assert os.path.isfile(os.path.join(mock_args.outdir, expected))


if __name__ == "__main__":
    pytest.main([__file__])
