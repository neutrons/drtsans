import pytest
from unittest.mock import Mock, patch
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, mtd
from drtsans.tof.eqsans.blocked_beam import subtract_blocked_beam


@pytest.fixture
def sample_workspace():
    """Create a sample input workspace"""
    CreateWorkspace(DataX=[1, 2, 3, 4, 5], DataY=[10, 20, 30, 40], DataE=[1, 1, 1, 1], OutputWorkspace="test_input_ws")
    yield "test_input_ws"
    if "test_input_ws" in mtd:
        DeleteWorkspace("test_input_ws")


@pytest.fixture
def blocked_beam_workspace():
    """Create a blocked beam raw workspace"""
    CreateWorkspace(
        DataX=[1, 2, 3, 4, 5], DataY=[2, 3, 4, 5], DataE=[0.5, 0.5, 0.5, 0.5], OutputWorkspace="blocked_beam_raw_histo"
    )
    yield "blocked_beam_raw_histo"
    for ws_name in ["blocked_beam_raw_histo", "blocked_beam_processed_histo"]:
        if ws_name in mtd:
            DeleteWorkspace(ws_name)


@pytest.fixture
def blocked_beam_mock():
    """Create a mock blocked_beam object"""
    blocked_beam = Mock()
    blocked_beam.data = "blocked_beam_raw_histo"
    return blocked_beam


@patch("drtsans.tof.eqsans.blocked_beam.subtract_dark_current")
@patch("drtsans.tof.eqsans.blocked_beam.normalize_by_flux")
@patch("drtsans.tof.eqsans.blocked_beam.subtract_background")
def test_processes_blocked_beam_workspace_first_time(
    mock_subtract_bg, mock_normalize, mock_subtract_dc, sample_workspace, blocked_beam_workspace, blocked_beam_mock
):
    """Test that blocked beam workspace is cloned and processed when it doesn't exist"""
    flux = Mock()
    subtract_blocked_beam(
        input_workspace=sample_workspace, blocked_beam=blocked_beam_mock, flux_method="proton charge", flux=flux
    )

    assert "blocked_beam_processed_histo" in mtd
    mock_normalize.assert_called_once_with("blocked_beam_processed_histo", flux, method="proton charge")
    mock_subtract_dc.assert_not_called()
    mock_subtract_bg.assert_called_once()


@patch("drtsans.tof.eqsans.blocked_beam.normalize_by_flux")
@patch("drtsans.tof.eqsans.blocked_beam.subtract_background")
def test_reuses_existing_processed_workspace(
    mock_subtract_bg, mock_normalize, sample_workspace, blocked_beam_workspace, blocked_beam_mock
):
    """Test that existing processed workspace is reused without reprocessing"""
    CreateWorkspace(
        DataX=[1, 2, 3, 4, 5],
        DataY=[1, 1, 1, 1],
        DataE=[0.1, 0.1, 0.1, 0.1],
        OutputWorkspace="blocked_beam_processed_histo",
    )
    flux = Mock()
    subtract_blocked_beam(
        input_workspace=sample_workspace, blocked_beam=blocked_beam_mock, flux_method="proton charge", flux=flux
    )

    mock_normalize.assert_not_called()
    mock_subtract_bg.assert_called_once()


@patch("drtsans.tof.eqsans.blocked_beam.subtract_dark_current")
@patch("drtsans.tof.eqsans.blocked_beam.normalize_by_flux")
@patch("drtsans.tof.eqsans.blocked_beam.subtract_background")
def test_processes_with_dark_current(
    mock_subtract_bg, mock_normalize, mock_subtract_dc, sample_workspace, blocked_beam_workspace, blocked_beam_mock
):
    """Test that dark current is subtracted when provided"""
    dark_current = Mock()
    dark_current.data = "dark_current_ws"
    flux = Mock()

    subtract_blocked_beam(
        input_workspace=sample_workspace,
        blocked_beam=blocked_beam_mock,
        flux_method="proton charge",
        flux=flux,
        dark_current=dark_current,
    )
    mock_subtract_dc.assert_called_once_with("blocked_beam_processed_histo", "dark_current_ws")

    mock_normalize.assert_called_once()
    mock_subtract_bg.assert_called_once()
