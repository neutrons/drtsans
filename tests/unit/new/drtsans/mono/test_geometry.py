import pytest
from os.path import join as path_join
from mantid.simpleapi import LoadEventNexus
from drtsans.settings import unique_workspace_dundername
from drtsans.mono.geometry import beam_radius


def test_beam_radius(reference_dir):
    """Test beam radius calculation for mono SANS instruments (BIO and GP-SANS)

    Parameters
    ----------
    reference_dir: str
        directory for test data

    Returns
    -------

    """
    workspace = LoadEventNexus(Filename=path_join(reference_dir.new.gpsans, 'geometry', 'CG2_7614.nxs.h5'),
                               OutputWorkspace=unique_workspace_dundername(), MetaDataOnly=True, LoadLogs=True)
    assert beam_radius(workspace, unit='mm') == pytest.approx(16.0, abs=0.1)
    workspace.delete()
