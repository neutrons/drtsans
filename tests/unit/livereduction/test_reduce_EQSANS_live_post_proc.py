import io
import logging
import pathlib
import shutil
import sys
import types
from unittest import mock

from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import DeleteWorkspace, MoveInstrumentComponent, mtd, Rebin
import numpy as np
import pytest

from drtsans.instruments import empty_instrument_workspace
from drtsans.path import load_module
from drtsans.samplelogs import SampleLogs
from drtsans.simulated_events import insert_background

_root_dir = pathlib.Path(__file__).parent.parent.parent.parent  # Go up 4 levels from test file


@pytest.fixture(scope="module")
def simulated_events() -> EventWorkspace:
    """Create a simulated EQSANS instrument workspace with 9 events per pixel

    Also adds sample logs:
    - start_time = 2023-08-01 00:00:00
    - run_number = 12345
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


def test_configure_error_buffer():
    """Test that configure_error_buffer creates a context manager that properly captures errors."""
    root = logging.getLogger()
    livescript = load_module(_root_dir / "scripts/livereduction/eqsans/reduce_EQSANS_live_post_proc.py")

    with livescript.configure_error_buffer() as buf:
        assert isinstance(buf, io.StringIO)
        # Verify handler was added
        handler = next(
            (h for h in root.handlers if isinstance(h, logging.StreamHandler) and getattr(h, "stream", None) is buf),
            None,
        )
        assert handler is not None
        assert handler.level == logging.ERROR
        assert handler.formatter is not None
        assert handler.formatter._fmt == "%(asctime)s - %(levelname)s - %(message)s"
        # Log a message at level below ERROR and verify it's not captured
        root.info("Test info message")
        captured = buf.getvalue()
        assert "Test info message" not in captured
        # Log an error and verify it's captured
        root.error("Test error message")
        captured = buf.getvalue()
        assert "Test error message" in captured
    # Verify cleanup: handler should be removed after exiting context
    handler_exists = any(
        isinstance(h, logging.StreamHandler) and getattr(h, "stream", None) is buf for h in root.handlers
    )
    assert not handler_exists, "Handler should be removed after context exit"


def test_livereduce(simulated_events, tmp_path):
    # mock the call to shutil.copytree inside livereduce to copy to tmp_path instead
    def mock_copytree_impl(src, dst, **kwargs):
        return shutil.copytree(src, str(tmp_path), **kwargs)

    mock_copytree = mock.MagicMock(side_effect=mock_copytree_impl)

    # the mock should behave as a context manager
    mock_add_to_sys_path = mock.MagicMock()
    mock_add_to_sys_path.__enter__ = lambda self, *_: self
    mock_add_to_sys_path.__exit__ = lambda *_: None

    mock_makedirs = mock.MagicMock(return_value=None)

    def mock_save_report(report, filename, logger):
        with open(filename, "w") as f:
            f.write(report)

    # Mock the reduce_EQSANS module functions
    mock_reduce_module = types.ModuleType("reduce_EQSANS")
    mock_reduce_module.reduce_events = mock.MagicMock(return_value="<report>")
    mock_reduce_module.footer = mock.MagicMock(return_value="<footer>")
    mock_reduce_module.save_report = mock_save_report
    mock_reduce_module.upload_report = mock.MagicMock()
    mock_reduce_module.LogContext = mock.MagicMock()
    sys.modules["reduce_EQSANS"] = mock_reduce_module

    # Mock SaveNexusProcessed to track calls
    mock_save_nexus = mock.MagicMock()

    try:
        # import there so mocking the module's symbols won't affect other tests
        livescript = load_module(_root_dir / "scripts/livereduction/eqsans/reduce_EQSANS_live_post_proc.py")
        livescript.add_to_sys_path = mock_add_to_sys_path
        livescript.events_file_exists = mock.MagicMock(return_value=False)
        livescript.copytree = mock_copytree
        livescript.makedirs = mock_makedirs
        livescript.SaveNexusProcessed = mock_save_nexus
        livescript.livereduce(simulated_events)

        # Verify SaveNexusProcessed was called with the events workspace
        assert mock_save_nexus.called, "SaveNexusProcessed should be called for live reduction"
        save_call_args = mock_save_nexus.call_args
        assert save_call_args[1]["InputWorkspace"] == simulated_events
        # Verify the temp file path contains the run number
        assert "EQSANS_12345_live.nxs" in save_call_args[1]["Filename"]

        # Verify reduce_events was called with sample_file parameter
        reduce_events_call = mock_reduce_module.reduce_events
        assert reduce_events_call.called
        call_kwargs = reduce_events_call.call_args[1]
        assert "sample_file" in call_kwargs
        assert "EQSANS_12345_live.nxs" in call_kwargs["sample_file"]

        assert (tmp_path / "EQSANS_12345.html").read_text() == "<report><footer>"
    finally:  # Clean up
        sys.modules.pop("reduce_EQSANS", None)


def test_livereduce_temp_file_created_in_temp_dir(simulated_events, tmp_path):
    """Test that the temporary sample file is created within the temporary directory."""
    temp_files_created = []

    def mock_save_nexus(InputWorkspace, Filename):
        temp_files_created.append(Filename)

    mock_add_to_sys_path = mock.MagicMock()
    mock_add_to_sys_path.__enter__ = lambda self, *_: self
    mock_add_to_sys_path.__exit__ = lambda *_: None

    mock_reduce_module = types.ModuleType("reduce_EQSANS")
    mock_reduce_module.reduce_events = mock.MagicMock(return_value="<report>")
    mock_reduce_module.footer = mock.MagicMock(return_value="<footer>")
    mock_reduce_module.save_report = mock.MagicMock()
    mock_reduce_module.upload_report = mock.MagicMock()
    mock_reduce_module.LogContext = mock.MagicMock()
    sys.modules["reduce_EQSANS"] = mock_reduce_module

    try:
        livescript = load_module(_root_dir / "scripts/livereduction/eqsans/reduce_EQSANS_live_post_proc.py")
        livescript.add_to_sys_path = mock_add_to_sys_path
        livescript.events_file_exists = mock.MagicMock(return_value=False)
        livescript.copytree = mock.MagicMock()
        livescript.makedirs = mock.MagicMock()
        livescript.SaveNexusProcessed = mock_save_nexus
        livescript.livereduce(simulated_events, publish=False)

        # Verify a temp file was created
        assert len(temp_files_created) == 1
        temp_file_path = temp_files_created[0]

        # Verify the temp file is in a temp directory (contains 'livereduce_' prefix)
        assert "livereduce_" in temp_file_path
        assert "EQSANS_12345_live.nxs" in temp_file_path
    finally:
        sys.modules.pop("reduce_EQSANS", None)


if __name__ == "__main__":
    pytest.main([__file__])
