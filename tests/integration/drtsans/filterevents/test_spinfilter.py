import os

from mantid.simpleapi import DeleteWorkspace, LoadEventNexus, mtd
import pytest

from drtsans.filterevents.spinfilter import SpinFilter
from drtsans.polarization import (
    SimulatedPolarizationLogs,
    TimesGeneratorSpecs,
)


class TestSpinFilter:
    @pytest.fixture(scope="function")
    def gpsans_workspace(self, datarepo_dir):
        """GPSANS run that is 300s long"""
        nexus_name = os.path.join(datarepo_dir.gpsans, "CG2_9166.nxs.h5")
        workspace = LoadEventNexus(
            Filename=nexus_name,
            LoadMonitors=False,
            OutputWorkspace=mtd.unique_hidden_name(),
        )
        return workspace

    @pytest.mark.datarepo
    def test_polarizer_flipper_analyzer_flipper(self, gpsans_workspace):
        """Test splitting events based on both polarizer and analyzer states, with flipper and veto.

        The analyzer interval (61 s) is chosen to be coprime with the polarizer interval (60 s)
        so that their state-change timestamps never coincide, avoiding the deduplication in
        SpinFilter._build_change_list that would otherwise silently drop all analyzer events.
        """
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "alive_duration": 1.0}),
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 61.0}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 61.0, "alive_duration": 1.0}),
        )
        logs.inject(workspace)
        sf = SpinFilter(workspace)
        sf.apply_filter("workspace_split")
        ws_group = mtd["workspace_split"]
        assert len(ws_group) == 4
        assert ws_group[0].getNumberEvents() == 883574
        assert ws_group[1].getNumberEvents() == 30957
        assert ws_group[2].getNumberEvents() == 863634
        assert ws_group[3].getNumberEvents() == 15511
        if mtd.doesExist("workspace_split"):
            DeleteWorkspace("workspace_split")

    @pytest.mark.datarepo
    def test_polarizer(self, gpsans_workspace):
        """Test splitting events based on polarizer state only."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "alive_duration": 1.0}),
        )
        logs.inject(workspace)
        sf = SpinFilter(workspace)
        sf.apply_filter("workspace_split")
        ws_group = mtd["workspace_split"]
        assert len(ws_group) == 2
        assert ws_group[0].getNumberEvents() == 1368108
        assert ws_group[1].getNumberEvents() == 910749
        if mtd.doesExist("workspace_split"):
            DeleteWorkspace("workspace_split")

    @pytest.mark.datarepo
    def test_analyzer(self, gpsans_workspace):
        """Test splitting events based on analyzer state only."""
        workspace = gpsans_workspace
        logs = SimulatedPolarizationLogs(
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 120}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 120.0, "alive_duration": 2.0}),
        )
        logs.inject(workspace)
        sf = SpinFilter(workspace)
        sf.apply_filter("workspace_split")
        ws_group = mtd["workspace_split"]
        assert len(ws_group) == 2
        assert ws_group[0].getNumberEvents() == 1373190
        assert ws_group[1].getNumberEvents() == 909694
        if mtd.doesExist("workspace_split"):
            DeleteWorkspace("workspace_split")

    @pytest.mark.datarepo
    def test_no_filtering(self, gpsans_workspace):
        """
        Verify that if no polarization devices are present SpinFilter raises ValueError.

        The new SpinFilter.__init__ eagerly checks for active polarizer/analyzer logs and
        raises ValueError immediately when neither is found, rather than deferring to
        apply_filter and returning the raw workspace unchanged.
        """
        workspace = gpsans_workspace
        with pytest.raises(ValueError, match="No active polarizer or analyzer"):
            SpinFilter(workspace)
