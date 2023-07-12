# local imports
from drtsans.mono.transmission import calculate_transmission
from drtsans.mono.biosans import prepare_data
from drtsans.settings import unique_workspace_dundername

# third party imports
from mantid.simpleapi import DeleteWorkspaces, LoadNexusProcessed
from numpy.testing import assert_almost_equal
import pytest

# standard imports
from unittest.mock import patch as mock_patch
import os


def _mock_LoadEventNexus(*args, **kwargs):
    # Substitute LoadEventNexus with LoadNexusProcessed because our synthetic files were created with SaveNexus
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


@pytest.mark.skip(reason="Skip until the mantid nightly conda package is updated, with timestamp > 2030-07-12")
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
@mock_patch("drtsans.load.__monitor_counts")
def test_transmission(mock_monitor_counts, biosans_synthetic_dataset):
    data = biosans_synthetic_dataset
    mock_monitor_counts.return_value = data["monitor_counts"]
    empty, sample = [
        prepare_data(
            data=os.path.join(biosans_synthetic_dataset["data_dir"], f"CG3_{data['runs'][run]}.nxs.h5"),
            mask_detector="wing_detector, midrange_detector",
            btp={"Pixel": "1-18,239-256"},
            centering_options={"IntegrationRadius": 0.03},
            detector_offset=0,
            sample_offset=0,
            center_x=0,
            center_y=0,
            flux_method="monitor",
            solid_angle=False,
            sample_thickness=0.1,
            output_workspace=unique_workspace_dundername(),
        )
        for run in ("empty_transmission", "sample_transmission")
    ]
    transmission = calculate_transmission(sample, empty)
    assert_almost_equal(transmission.extractY()[0][0], 0.75, decimal=2)
    DeleteWorkspaces([sample, empty, transmission])


if __name__ == "__main__":
    pytest.main([__file__])
