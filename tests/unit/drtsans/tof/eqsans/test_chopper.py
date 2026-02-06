import pytest
from mantid.simpleapi import AddSampleLog, CreateSingleValuedWorkspace, mtd
from os.path import join as pjn
from numpy.testing import assert_almost_equal

from drtsans.tof.eqsans.chopper import (
    EQSANSFourChoppersConfiguration,
    EQSANSSixChoppersConfiguration,
    FrameMode,
    EQSANSDiskChopperSet,
)


class TestEQSANSDiskChopperSet:
    @pytest.mark.datarepo
    def test_transmitted_bands(self, datarepo_dir, clean_workspace):
        # Test transmitted bands in skipping mode
        file_name = pjn(datarepo_dir.eqsans, "test_chopper", "EQSANS_92353_no_events.nxs")
        chs = EQSANSDiskChopperSet(file_name)
        clean_workspace("ws")
        assert chs.frame_mode == FrameMode.skip
        # prompt pulse
        wb = chs.transmission_bands(pulsed=True)
        assert len(wb) == 1
        assert_almost_equal((wb[0].min, wb[0].max), (2.48, 6.13), decimal=2)
        # skipped pulse has a delay of chs.period/2 (16666.6 micro-sec)
        wb = chs.transmission_bands(delay=chs.period / 2, pulsed=True)
        assert len(wb) == 1
        assert_almost_equal((wb[0].min, wb[0].max), (9.64, 13.41), decimal=2)

        # Test transmitted bands in non-skipping mode
        #
        # porasil 1m
        file_name = pjn(datarepo_dir.eqsans, "test_chopper", "EQSANS_92164_no_events.nxs")
        chs = EQSANSDiskChopperSet(file_name)
        assert chs.frame_mode == FrameMode.not_skip
        wb = chs.transmission_bands(pulsed=True)
        assert len(wb) == 2  # there's a small leakage, the second band
        assert_almost_equal((wb[0].min, wb[0].max), (2.47, 6.66), decimal=2)
        # previous pulse has a delay of chs.period (16666.6 micro-sec)
        wb = chs.transmission_bands(delay=chs.period, pulsed=True)
        assert len(wb) == 1  # the small leakage from the previous pulse
        assert_almost_equal((wb[0].min, wb[0].max), (27.58, 27.60), decimal=2)
        #
        # porasil 4m
        file_name = pjn(datarepo_dir.eqsans, "test_chopper", "EQSANS_92149_no_events.nxs")
        chs = EQSANSDiskChopperSet(file_name)
        assert chs.frame_mode == FrameMode.not_skip
        wb = chs.transmission_bands(pulsed=True)
        assert len(wb) == 1
        assert_almost_equal((wb[0].min, wb[0].max), (9.91, 13.63), decimal=2)
        # skipped pulse has a delay of chs.period (16666.6 micro-sec)
        wb = chs.transmission_bands(delay=chs.period, pulsed=True)
        assert len(wb) == 1  # we are working in the second frame
        assert_almost_equal((wb[0].min, wb[0].max), (9.91, 13.63), decimal=2)
        #
        # porasil 8m
        file_name = pjn(datarepo_dir.eqsans, "test_chopper", "EQSANS_92144_no_events.nxs")
        chs = EQSANSDiskChopperSet(file_name)
        assert chs.frame_mode == FrameMode.not_skip
        wb = chs.transmission_bands(pulsed=True)
        assert len(wb) == 1
        assert_almost_equal((wb[0].min, wb[0].max), (11.89, 14.98), decimal=2)
        # skipped pulse has a delay of chs.period (16666.6 micro-sec)
        wb = chs.transmission_bands(delay=chs.period, pulsed=True)
        assert len(wb) == 1  # we are working in the second frame
        assert_almost_equal((wb[0].min, wb[0].max), (11.89, 14.98), decimal=2)

    def test_get_chopper_configuration(self):
        """Test getting the chopper configuration depending on the available sample logs."""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="Speed1", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase1", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Speed2", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase2", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Speed3", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase3", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Speed4", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase4", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="frequency", LogText="60", LogType="Number Series")

        # 4 choppers
        chs = EQSANSDiskChopperSet(workspace)
        assert chs._n_choppers == 4
        assert isinstance(chs.chopper_config, EQSANSFourChoppersConfiguration)

        AddSampleLog(Workspace=workspace, LogName="Speed5", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase5", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Speed6", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="Phase6", LogText="1", LogType="Number Series")

        # 6 choppers
        chs = EQSANSDiskChopperSet(workspace)
        assert chs._n_choppers == 6
        assert isinstance(chs.chopper_config, EQSANSSixChoppersConfiguration)


if __name__ == "__main__":
    pytest.main([__file__])
