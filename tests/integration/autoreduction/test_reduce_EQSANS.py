import pathlib
import os
from unittest.mock import Mock

import pytest

from drtsans.path import load_module

# Add the root directory to the path
_root_dir = pathlib.Path(__file__).parent.parent.parent.parent  # Go up 4 levels from test file
reduce_EQSANS = load_module(_root_dir / "scripts/autoreduction/reduce_EQSANS.py")


@pytest.mark.mount_eqsans
def test_autoreduce_sample(tmp_path):
    """The counts for this run are not isotropic, so we can use it to assert the I(Qx, Qy)
    plot reflects the asymmetry. Here's a rough ASCII representation of what the
    I(Qx, Qy) plot should look like. The center of the intensity maximum is shifted
    towards negative Qy values.
    Qy
    ^
    |..........................................
    |..........................................
    |..........................................
    |....................:::::.................
    --------------------:::::::------------------> Qx
    |..................:::::::::...............
    |..................:::::::::...............
    |....................:::::.................
    |..........................................
    """
    mock_args = Mock()
    mock_args.events_file = "/SNS/EQSANS/IPTS-35884/nexus/EQSANS_172835.nxs.h5"
    mock_args.outdir = str(tmp_path)
    mock_args.no_publish = True  # Disable publishing to the live data server
    reduce_EQSANS.autoreduce(mock_args)
    filenames = [
        "autoreduce_172835.log",
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
    for expected in filenames:
        assert os.path.isfile(os.path.join(mock_args.outdir, expected))


if __name__ == "__main__":
    pytest.main([__file__])
