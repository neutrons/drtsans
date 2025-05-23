import os
from unittest import mock

import pytest
from mantid.simpleapi import mtd, LoadEventNexus, DeleteWorkspace

from drtsans.filter_events import split_events
from drtsans.polarization import (
    SimulatedPolarizationLogs,
    TimesGeneratorSpecs,
)


class TestFilterEvents:
    @pytest.fixture(scope="function")
    def gpsans_workspace(self, datarepo_dir):
        """GP-SANS run that is 300 s long"""
        nexus_name = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
        workspace = LoadEventNexus(
            Filename=nexus_name,
            LoadMonitors=False,
            OutputWorkspace=mtd.unique_hidden_name(),
        )
        return workspace

    @pytest.mark.datarepo
    def test_polarizer_flipper_analyzer_flipper(self, gpsans_workspace):
        """Test splitting events based on both polarizer and analyzer states, with flipper and veto."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "veto_duration": 1.0}),
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 120}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 120.0, "veto_duration": 2.0}),
        )
        logs.inject(workspace)
        ws_group = split_events("workspace_split", input_workspace=workspace)
        assert len(ws_group) == 4
        assert ws_group[0].getNumberEvents() == 912779
        assert ws_group[1].getNumberEvents() == 452402
        assert ws_group[2].getNumberEvents() == 455174
        assert ws_group[3].getNumberEvents() == 450732
        if mtd.doesExist(str(ws_group)):
            DeleteWorkspace(ws_group)

    @pytest.mark.datarepo
    def test_polarizer(self, gpsans_workspace):
        """Test splitting events based on polarizer state."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "veto_duration": 1.0}),
        )
        logs.inject(workspace)
        ws_group = split_events("workspace_split", input_workspace=workspace)
        assert len(ws_group) == 2
        assert ws_group[0].getNumberEvents() == 1367953
        assert ws_group[1].getNumberEvents() == 910789
        if mtd.doesExist(str(ws_group)):
            DeleteWorkspace(ws_group)

    @pytest.mark.datarepo
    def test_analyzer(self, gpsans_workspace):
        """Test splitting events based on analyzer state."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 120}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 120.0, "veto_duration": 2.0}),
        )
        logs.inject(workspace)
        ws_group = split_events("workspace_split", input_workspace=workspace)
        assert len(ws_group) == 2
        assert ws_group[0].getNumberEvents() == 1373177
        assert ws_group[1].getNumberEvents() == 909718
        if mtd.doesExist(str(ws_group)):
            DeleteWorkspace(ws_group)

    @pytest.mark.datarepo
    def test_split_with_input_file(self, gpsans_workspace):
        """Test splitting events using an input Nexus file."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "veto_duration": 1.0}),
        )
        logs.inject(workspace)
        with mock.patch("drtsans.filter_events.LoadEventNexus") as mock_loadeventnexus:
            mock_loadeventnexus.return_value = workspace
            ws_group = split_events("workspace_split", file_path="/path/to/file")
        assert len(ws_group) == 2
        assert ws_group[0].getNumberEvents() == 1367953
        assert ws_group[1].getNumberEvents() == 910789
        if mtd.doesExist(str(ws_group)):
            DeleteWorkspace(ws_group)

    @pytest.mark.datarepo
    def test_no_filtering(self, gpsans_workspace):
        """
        Verify that if no filtering is requested, the original workspace is returned unchanged.
        """
        workspace = gpsans_workspace
        num_events = workspace.getNumberEvents()
        ws_group = split_events("workspace_split", input_workspace=workspace)
        assert len(ws_group) == 1
        assert ws_group[0].getNumberEvents() == num_events
        if mtd.doesExist(str(ws_group)):
            DeleteWorkspace(ws_group)
