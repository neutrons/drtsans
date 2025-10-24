import os
from unittest.mock import Mock

import pytest

from .script_locator import reduce_EQSANS


@pytest.mark.mount_eqsans
def test_autoreduce_sample(simulated_events, tmp_path):
    mock_args = Mock()
    mock_args.events_file = "/SNS/EQSANS/IPTS-34577/nexus/EQSANS_162568.nxs.h5"
    mock_args.outdir = str(tmp_path)
    mock_args.no_publish = True  # Disable publishing to the live data server
    reduce_EQSANS.autoreduce(mock_args)
    filenames = [
        "autoreduce_162568.log",
        "EQSANS_162568.html",
        "EQSANS_162568_Iq.dat",
        "EQSANS_162568_Iq.png",
        "EQSANS_162568_Iqxqy.dat",
        "EQSANS_162568_Iqxqy.h5",
        "EQSANS_162568_Iqxqy.png",
        "EQSANS_162568_processed.nxs",
        "EQSANS_162568_reduction_log.hdf",
        "reduction_options_162568.json",
    ]
    for expected in filenames:
        assert os.path.isfile(os.path.join(mock_args.outdir, expected))


if __name__ == "__main__":
    pytest.main([__file__])
