import pytest
from os.path import join as path_join
from mantid.simpleapi import LoadEventNexus, mtd
from drtsans.mono.geometry import beam_radius


@pytest.mark.datarepo
def test_beam_radius(datarepo_dir):
    """Test beam radius calculation for mono SANS instruments (BIO and GP-SANS)"""

    workspace = LoadEventNexus(
        Filename=path_join(datarepo_dir.gpsans, "CG2_9188.nxs.h5"),
        OutputWorkspace=mtd.unique_hidden_name(),
        MetaDataOnly=True,
        LoadLogs=True,
    )
    assert beam_radius(workspace, unit="mm") == pytest.approx(34.1, abs=0.1)  # calculated value = 34.137181990397
    workspace.delete()


if __name__ == "__main__":
    pytest.main([__file__])
