import pytest
import numpy as np

from drtsans.mono.dark_current import normalize_dark_current
from drtsans.samplelogs import SampleLogs


# Equations from the test, thus I_dc is the array of dark current intensities which we must compare to
x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
I_dc = (x + y) * 2. + 5.
t_sam = 5.  # sample measure for 100s
t_dc = 3600.  # dark current collected for 1 hr.


@pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 5, 'Ny': 5}], indirect=True)
def test_dark_current_normalize_to_workspace(workspace_with_instrument):
    """Test of dark current normalization
    For details see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/156
    and also https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/174

    dev - Jose Borreguero <borreguerojm@ornl.gov>, Steven Hahn <hahnse@ornl.gov>, Jiao Lin <linjiao@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>

    **drtsans functions used:**
    ~drtsans.samplelogs.SampleLogs
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
    ~drtsans.mono.dark_current.normalize_dark_current
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/dark_current.py>
    """
    # Create sample workspace with irrelevant intensities. We need to insert the duration of the sample run as one of
    # the log entries in the sample workspace. We have also associated a detector panel, more realistic albeit
    # unnecessary for the purposes of normalizing the dark current
    data_workspace = workspace_with_instrument(intensities=np.zeros((5, 5)), view='pixel')
    SampleLogs(data_workspace).insert('duration', t_sam, 'second')

    # Create dark current workspace, insert the duration of the dark current run as one of the log entries in the
    # dark current workspace.
    dark_current_workspace = workspace_with_instrument(intensities=np.array(I_dc), view='pixel')
    SampleLogs(dark_current_workspace).insert('duration', t_dc, 'second')

    # Call our normalization function on the dark current workspace. Notice that normalization by the duration of
    # the sample duration is not carried out. The reason is that we can re-use this normalized dark current with
    # sample runs of different duration.
    normalize_dark_current(dark_current_workspace)

    I_dcnorm_step = t_sam / t_dc * I_dc  # normalized current following the equations of the test

    # Compare with our normalized dark current. Notice that at this step we use the duration of the sample run
    assert np.allclose(t_sam * dark_current_workspace.extractY().ravel(), I_dcnorm_step.ravel())


if __name__ == '__main__':
    pytest.main([__file__])
