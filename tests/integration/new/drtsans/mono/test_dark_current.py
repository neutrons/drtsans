import pytest
import numpy as np

from drtsans.mono.dark_current import normalize_dark_current
from drtsans.samplelogs import SampleLogs

x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
I_dc = (x + y) * 2. + 5.
t_sam = 5.  # sample measure for 100s
t_dc = 3600.  # dark current collected for 1 hr.


@pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 5, 'Ny': 5}], indirect=True)
def test_dark_current_normalize_to_workspace(workspace_with_instrument):
    data_ws = workspace_with_instrument(intensities=np.zeros((5, 5)), view='pixel')
    SampleLogs(data_ws).insert('duration', t_sam, 'second')
    dark_ws = workspace_with_instrument(intensities=np.array(I_dc),
                                        view='pixel')
    SampleLogs(dark_ws).insert('duration', t_dc, 'second')
    normalize_dark_current(dark_ws)
    I_dcnorm_step = t_sam / t_dc * I_dc
    assert np.allclose(dark_ws.extractY().ravel(), I_dcnorm_step.ravel())


if __name__ == '__main__':
    pytest.main()
