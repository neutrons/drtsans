import pathlib
from unittest import mock
import shutil

from mantid.simpleapi import LoadEventNexus, mtd
import pytest

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


if __name__ == "__main__":
    pytest.main([__file__])
