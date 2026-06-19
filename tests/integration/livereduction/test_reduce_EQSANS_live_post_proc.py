import pathlib
from unittest import mock
import shutil

from mantid.simpleapi import LoadEventNexus, SaveNexusProcessed, mtd
import pytest

from drtsans.load import load_events
from drtsans.path import load_module

# Add the repo's root directory to the path
_root_dir = pathlib.Path(__file__).parent.parent.parent.parent  # Go up 4 levels from test file
livescript = load_module(_root_dir / "scripts/livereduction/eqsans/reduce_EQSANS_live_post_proc.py")
reduce_EQSANS_posixpath = _root_dir / "scripts/autoreduction/reduce_EQSANS.py"


@pytest.mark.mount_eqsans
def test_livereduce(tmp_path):
    events = LoadEventNexus(
        Filename="/SNS/EQSANS/IPTS-35884/nexus/EQSANS_172835.nxs.h5", OutputWorkspace=mtd.unique_hidden_name()
    )
    livescript.GLOBAL_AR_DIR = str(reduce_EQSANS_posixpath.parent)  # reduce_EQSANS.py in this codebase
    # mock the calls to os.makedirs and shutil.copytree inside livereduce
    livescript.makedirs = mock.MagicMock(return_value=None)

    def mock_copytree_impl(src, dst, **kwargs):
        return shutil.copytree(src, str(tmp_path), **kwargs)

    livescript.copytree = mock.MagicMock(side_effect=mock_copytree_impl)
    livescript.events_file_exists = mock.MagicMock(return_value=False)
    livescript.livereduce(events, publish=False)
    expected_files = [
        "EQSANS_172835.html",
        "EQSANS_172835_Iq.dat",
        "EQSANS_172835_Iq.png",
        "EQSANS_172835_Iqxqy.dat",
        "EQSANS_172835_Iqxqy.h5",
        "EQSANS_172835_Iqxqy.png",
        "EQSANS_172835_processed.nxs",
        "EQSANS_172835_reduction_log.hdf",
        "reduction_options_172835.json",
    ]

    for expected in expected_files:
        assert (tmp_path / expected).is_file(), f"{expected} was not created."


@pytest.mark.mount_eqsans
def test_save_and_load_processed_nexus_round_trip(tmp_path):
    """Test that SaveNexusProcessed + load_events with allow_processed_nexus works.

    This test verifies the round-trip of saving an EventWorkspace to a processed
    Nexus file and loading it back using the new allow_processed_nexus functionality.
    This is the core mechanism used for live reduction.
    """
    # Load original events
    original_ws_name = mtd.unique_hidden_name()
    original_events = LoadEventNexus(
        Filename="/SNS/EQSANS/IPTS-35884/nexus/EQSANS_172835.nxs.h5", OutputWorkspace=original_ws_name
    )

    # Save to a processed Nexus file (simulating what live reduction does)
    temp_file = str(tmp_path / "EQSANS_172835_live.nxs")
    SaveNexusProcessed(InputWorkspace=original_events, Filename=temp_file)

    # Load back using the new allow_processed_nexus parameter
    loaded_ws_name = mtd.unique_hidden_name()
    loaded_events = load_events(
        run=temp_file,
        output_workspace=loaded_ws_name,
        allow_processed_nexus=True,
    )

    # Verify the loaded workspace has the same number of histograms
    assert loaded_events.getNumberHistograms() == original_events.getNumberHistograms()

    # Clean up
    mtd.remove(original_ws_name)
    mtd.remove(loaded_ws_name)


if __name__ == "__main__":
    pytest.main([__file__])
