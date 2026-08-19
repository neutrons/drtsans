from copy import deepcopy
import json
import os
import pathlib
import tempfile
from unittest.mock import patch, MagicMock

from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import DeleteWorkspace, MoveInstrumentComponent, mtd, Rebin
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

from drtsans.instruments import empty_instrument_workspace
from drtsans.path import load_module
from drtsans.samplelogs import SampleLogs
from drtsans.simulated_events import insert_background

_root_dir = pathlib.Path(__file__).parent.parent.parent.parent  # Go up 4 levels from test file
reduce_EQSANS = load_module(_root_dir / "scripts" / "autoreduction" / "reduce_EQSANS.py")


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
    rng = np.random.default_rng(7495230093183)  # add seed for deterministic results in tests
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
        lambda_distribution=lambda n_events: rng.normal(loc=5.0, scale=0.1, size=n_events),
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
    assert reduce_EQSANS.CONDA_ENV in ("sans", "sans_qa", "sans_dev")


def test_upload_report():
    """Test successful plot upload"""
    mock_logger = MagicMock()
    with patch.object(reduce_EQSANS, "publish_plot") as mock_publish_plot:
        mock_publish_plot.return_value = None
        reduce_EQSANS.upload_report("12345", "<div>test plot</div>", mock_logger)
        mock_publish_plot.assert_called_once_with("EQSANS", "12345", files={"file": "<div>test plot</div>"})


def test_save_report():
    """Test saving plot to file"""
    plot_div = "<div>test plot content</div>"
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".html") as temp_file:
        temp_filename = temp_file.name
    try:
        mock_logger = MagicMock()
        reduce_EQSANS.save_report(plot_div, temp_filename, mock_logger)
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
        mock_logger = MagicMock()
        reduce_EQSANS.save_report(plot_div, temp_filename, mock_logger)
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


def test_autoreduce_rejects_empty_event_workspace(tmp_path, mocker):
    events_file = tmp_path / "EQSANS_186603.nxs.h5"
    events_file.touch()
    empty_events = MagicMock()
    empty_events.getNumberEvents.return_value = 0
    mocker.patch.object(reduce_EQSANS, "LoadEventNexus", return_value=empty_events)

    args = MagicMock(events_file=str(events_file), outdir=str(tmp_path), no_publish=True)
    with pytest.raises(RuntimeError, match="has no events after loading with LoadEventNexus"):
        reduce_EQSANS.autoreduce(args)


def test_intensity_array(simulated_events):
    x, y, z = reduce_EQSANS.intensity_array(simulated_events)
    assert z.shape == (reduce_EQSANS.PIXELS_PER_TUBE, reduce_EQSANS.TUBES_IN_DETECTOR1)
    # event count in the first pixel, but all pixels should have the same count
    count = np.sum(simulated_events.readY(0))
    assert_almost_equal(np.average(z.data[~z.mask]), np.log(count), decimal=3)


@pytest.mark.parametrize("options_location", ["output_dir", "ipts", "shared"])
def test_reduce_sample_forces_fallback_beam_center(simulated_events, tmp_path, monkeypatch, mocker, options_location):
    """An unattended reduction turns the fallback beam center on, whichever input JSON it settles on

    The input file below asks for `false`, so a passing test also shows the amendment overrides it.
    """
    monkeypatch.chdir(tmp_path)  # the script chdir's to AUTOREDUCE_DIR; restore the cwd afterwards
    run_number = str(simulated_events.getRunNumber())
    output_dir, shared_dir, ipts_dir = tmp_path / "output", tmp_path / "shared", tmp_path / "IPTS-12345"
    for directory in (output_dir, shared_dir, ipts_dir):
        directory.mkdir()
    monkeypatch.setattr(reduce_EQSANS, "AUTOREDUCE_DIR", str(shared_dir))
    monkeypatch.setattr(reduce_EQSANS, "AUTOREDUCE_IPTS_DIR", str(tmp_path / "IPTS-{ipts}"))

    raw_options = {
        "instrumentName": "EQSANS",
        "iptsNumber": 12345,
        "sample": {"runNumber": run_number, "thickness": 1.0},
        "outputFileName": f"EQSANS_{run_number}",
        "configuration": {"outputDir": str(output_dir)},
        "beamCenter": {"runNumber": run_number, "useFallbackBeamCenter": False},
    }
    options_file = {
        "output_dir": output_dir / f"reduction_options_{run_number}.json",
        "ipts": ipts_dir / "reduction_options.json",
        "shared": shared_dir / "reduction_options.json",
    }[options_location]
    options_file.write_text(json.dumps(raw_options))

    # record the amendment, and skip validation because the dataSource validators need files on /SNS
    amendments = []
    real_update = reduce_EQSANS.update_reduction_parameters

    def record_amendment(parameters_original, parameter_changes, **kwargs):
        amendments.append(deepcopy(parameter_changes))
        return real_update(parameters_original, parameter_changes, validate=False, permissible=True)

    mocker.patch.object(reduce_EQSANS, "update_reduction_parameters", side_effect=record_amendment)
    mocker.patch.object(reduce_EQSANS, "load_all_files")
    mocker.patch.object(reduce_EQSANS, "reduce_single_configuration")
    mocker.patch.object(reduce_EQSANS, "plot_reduction_output")
    mocker.patch.object(reduce_EQSANS, "plotly_reduction_output", return_value="")
    mocker.patch.object(reduce_EQSANS, "reduce_non_sample", return_value="")
    mocker.patch.object(reduce_EQSANS, "GPR_AVAILABLE", False)

    reduce_EQSANS.reduce_sample(simulated_events, str(output_dir), MagicMock())

    assert amendments[0]["beamCenter"]["useFallbackBeamCenter"] is True
    if options_location == "shared":
        # the beam center run number joins the flag rather than replacing it
        assert amendments[0]["beamCenter"]["runNumber"] == run_number

    # the comprehensive options saved for the record carry the forced value, not the requested `false`
    saved = json.loads((output_dir / f"reduction_options_{run_number}.json").read_text())
    assert saved["beamCenter"]["useFallbackBeamCenter"] is True


if __name__ == "__main__":
    pytest.main([__file__])
