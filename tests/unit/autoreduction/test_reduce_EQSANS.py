import os
import requests
import tempfile
from unittest.mock import Mock, patch

from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import DeleteWorkspace, MoveInstrumentComponent, mtd, Rebin
import numpy as np
from numpy.testing import assert_almost_equal
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


def test_constants_values():
    """Test that constants have expected values"""
    assert reduce_EQSANS.TUBES_PER_EIGHTPACK == 8
    assert reduce_EQSANS.TUBES_IN_DETECTOR1 == 192
    assert reduce_EQSANS.PIXELS_PER_TUBE == 256
    assert reduce_EQSANS.CONDA_ENV == "sans"


def test_upload_report_success():
    """Test successful plot upload"""
    with patch.object(reduce_EQSANS, "publish_plot") as mock_publish_plot:
        mock_publish_plot.return_value = None
        reduce_EQSANS.upload_report("12345", "<div>test plot</div>")
        mock_publish_plot.assert_called_once_with("EQSANS", "12345", files={"file": "<div>test plot</div>"})


def test_upload_report_http_error():
    """Test handling of HTTP error during upload"""
    with patch.object(reduce_EQSANS, "publish_plot") as mock_publish_plot:
        mock_publish_plot.side_effect = requests.HTTPError("Test error")
        with patch.object(reduce_EQSANS, "logging") as mock_logging:
            reduce_EQSANS.upload_report("12345", "<div>test plot</div>")
            mock_logging.exception.assert_called_once_with("Publish plot failed with error Test error")


def test_save_report():
    """Test saving plot to file"""
    plot_div = "<div>test plot content</div>"
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".html") as temp_file:
        temp_filename = temp_file.name
    try:
        reduce_EQSANS.save_report(plot_div, temp_filename)
        with open(temp_filename, "r") as f:
            content = f.read()
        assert "<!DOCTYPE html>" in content
        assert plot_div in content
        assert "</html>" in content
    finally:
        os.unlink(temp_filename)


def test_save_report_with_html_structure():
    """Test that saved plot has correct HTML structure"""
    plot_div = "<div id='plotly-div'>Plot content</div>"
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".html") as temp_file:
        temp_filename = temp_file.name
    try:
        reduce_EQSANS.save_report(plot_div, temp_filename)
        with open(temp_filename, "r") as f:
            content = f.read()
        # Check HTML structure
        assert content.startswith("<!DOCTYPE html>")
        assert "<head>" in content
        assert "<title>Plotly Chart</title>" in content
        assert "<body>" in content
        assert content.strip().endswith("</html>")
    finally:
        os.unlink(temp_filename)


def test_parse_required_arguments():
    """Test parsing of required arguments"""
    test_args = ["test_events.nxs", "/output/dir"]

    with patch("sys.argv", ["script_name"] + test_args):
        args = reduce_EQSANS.parse_command_arguments()

    assert args.events_file == "test_events.nxs"
    assert args.outdir == "/output/dir"
    assert args.no_publish is False


def test_parse_all_arguments():
    """Test parsing of all arguments"""
    test_args = ["test_events.nxs", "/output/dir", "--no_publish"]

    with patch("sys.argv", ["script_name"] + test_args):
        args = reduce_EQSANS.parse_command_arguments()

    assert args.events_file == "test_events.nxs"
    assert args.outdir == "/output/dir"
    assert args.no_publish is True


def test_intensity_array(simulated_events):
    x, y, z = reduce_EQSANS.intensity_array(simulated_events)
    assert z.shape == (reduce_EQSANS.PIXELS_PER_TUBE, reduce_EQSANS.TUBES_IN_DETECTOR1)
    # event count in the first pixel, but all pixels should have the same count
    count = np.sum(simulated_events.readY(0))
    assert_almost_equal(np.average(z.data[~z.mask]), np.log(count), decimal=3)


@patch.object(reduce_EQSANS, "LoadEventNexus")
@patch.object(reduce_EQSANS, "parse_command_arguments")
@patch.object(reduce_EQSANS, "intensity_array")
@patch.object(reduce_EQSANS, "plot_heatmap")
@patch.object(reduce_EQSANS, "upload_report")
@patch.object(reduce_EQSANS, "save_report")
def test_autoreduce_with_publish_and_save(
    mock_save_report,
    mock_upload_report,
    mock_plot_heatmap,
    mock_intensity_array,
    mock_parse_args,
    mock_load,
    simulated_events,
    tmp_path,
):
    """Test autoreduce function with both publish and save enabled"""
    # Setup mocks
    mock_load.return_value = simulated_events

    mock_args = Mock()
    events_file = tmp_path / "invalid.nxs"
    events_file.touch()
    mock_args.events_file = str(events_file)
    mock_args.outdir = str(tmp_path)
    mock_args.no_publish = False
    mock_parse_args.return_value = mock_args

    mock_intensity_array.return_value = (
        np.array([1, 2, 3]),  # x values
        np.array([1, 2, 3]),  # y values
        np.ma.array([[1, 2], [3, 4]]),  # z values (masked array)
    )

    mock_plot_heatmap.return_value = "<div>plot content</div>"

    reduce_EQSANS.autoreduce()

    mock_intensity_array.assert_called_once()
    mock_plot_heatmap.assert_called_once()
    mock_upload_report.assert_called_once_with(12345, "<div>plot content</div>")
    mock_save_report.assert_called_once()


if __name__ == "__main__":
    pytest.main([__file__])
