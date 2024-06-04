# local imports
from drtsans.mono.biosans import find_beam_center, prepare_data

# third party imports
from mantid.simpleapi import DeleteWorkspace, LoadNexusProcessed, mtd
from numpy.testing import assert_allclose
import pytest

# standard imports
from unittest.mock import patch as mock_patch
import os


def _mock_LoadEventAsWorkspace2D(*args, **kwargs):
    # Substitute LoadEventAsWorkspace2D with LoadNexusProcessed because our synthetic files were created with SaveNexus
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


@pytest.mark.datarepo
@mock_patch("drtsans.load.LoadEventAsWorkspace2D", new=_mock_LoadEventAsWorkspace2D)
@mock_patch("drtsans.load.__monitor_counts")
def test_biosans_find_beam_center(mock_monitor_counts, biosans_synthetic_dataset):
    mock_monitor_counts.return_value = biosans_synthetic_dataset["monitor_counts"]
    beam_center_ws = prepare_data(
        data=os.path.join(biosans_synthetic_dataset["data_dir"], "CG3_92300.nxs.h5"),
        mask_detector="wing_detector, midrange_detector",
        btp={"Pixel": "1-18,239-256"},
        centering_options={"IntegrationRadius": 0.03},
        detector_offset=0,
        sample_offset=0,
        center_x=0,
        center_y=0,
        flux_method="monitor",
        solid_angle=False,
        output_workspace=mtd.unique_hidden_name(),
        sample_thickness=0.1,
    )
    center_x, center_y, center_y_wing, center_y_midrange, _ = find_beam_center(beam_center_ws)
    assert_allclose([center_x, center_y], [-0.008, -0.023], atol=0.004)
    DeleteWorkspace(beam_center_ws)


if __name__ == "__main__":
    pytest.main([__file__])
