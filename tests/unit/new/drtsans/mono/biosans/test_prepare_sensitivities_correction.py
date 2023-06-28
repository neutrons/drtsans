# local imports
from drtsans.mono.biosans.prepare_sensitivities_correction import PrepareSensitivityCorrection
from drtsans.path import abspath
from drtsans.settings import amend_config

# third party imports
import pytest
from mantid.simpleapi import LoadNexusProcessed

# standard imports
from unittest.mock import patch as mock_patch


class TestPrepareSensitivityCorrection:
    def test_init(self):
        preparer_main = PrepareSensitivityCorrection(component="detector1")
        assert preparer_main._instrument == "CG3"
        assert preparer_main._component == "detector1"
        preparer_midrange = PrepareSensitivityCorrection(component="midrange_detector")
        assert preparer_midrange._instrument == "CG3"
        assert preparer_midrange._component == "midrange_detector"

    def test_curved_detectors(self, biosans_synthetic_dataset):
        preparer_midrange = PrepareSensitivityCorrection(component="midrange_detector")
        preparer_midrange._flood_runs = [biosans_synthetic_dataset["flood"]]
        with amend_config(data_dir=str(biosans_synthetic_dataset["data_dir"])):
            assert preparer_midrange.curved_detectors == ["wing_detector", "midrange_detector"]

    def test_prepare_data_opts(self):
        preparer_midrange = PrepareSensitivityCorrection(component="midrange_detector")
        beam_center = [0.0, 1.0, 2.0, 3.0]
        opts_expected = dict(
            center_x=0.0, center_y=1.0, flux_method="monitor", overwrite_instrument=False, enforce_use_nexus_idf=False
        )
        opts = preparer_midrange._prepare_data_opts(beam_center[:2])
        assert opts == opts_expected

        opts_expected["center_y_wing"] = 2.0
        opts = preparer_midrange._prepare_data_opts(beam_center[:3])
        assert opts == opts_expected

        opts_expected["center_y_midrange"] = 3.0
        opts = preparer_midrange._prepare_data_opts(beam_center)
        assert opts == opts_expected

    @mock_patch("drtsans.load.LoadEventNexus")
    @mock_patch("drtsans.load.__monitor_counts")
    def test_get_beam_center_workspace(
        self,
        mock_monitor_counts,
        mock_load_LoadEventNexus,
        biosans_synthetic_dataset,
        clean_workspace,
    ):
        with amend_config(data_dir=str(biosans_synthetic_dataset["data_dir"])):
            # mock loading the beam center run and the monitor count
            beam_center_run = biosans_synthetic_dataset["beam_center"]
            workspace_beam_center = LoadNexusProcessed(
                Filename=abspath(beam_center_run, instrument="CG3"), OutputWorkspace=f"BC_CG3_{beam_center_run}"
            )
            clean_workspace(workspace_beam_center)  # mark for deletion upon exit or exception
            mock_load_LoadEventNexus.return_value = workspace_beam_center
            mock_monitor_counts.return_value = 42
            # test the method
            preparer_main = PrepareSensitivityCorrection(component="detector1")
            preparer_main._extra_mask_dict = dict(Components=["wing_detector", "midrange_detector"])
            preparer_main._flood_runs = [biosans_synthetic_dataset["flood"]]
            prepared_beam_center = preparer_main._get_beam_center_workspace(f"CG3_{beam_center_run}")
            clean_workspace(prepared_beam_center)  # mark for deletion upon exit or exception


if __name__ == "__main__":
    pytest.main([__file__])
