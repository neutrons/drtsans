import pytest
import numpy as np

x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
I_dc = (x + y) * 2. + 5.
t_sam = 5  # sample measure for 100s
t_dc = 3600  # dark current collected for 1 hr.
np.set_printoptions(precision=4, suppress=True)


@pytest.mark.parametrize('generic_workspace', [{'Nx': 5, 'Ny': 5}], indirect=True)
def test_dark_current_subtract_dark_current(generic_workspace):
    ws = generic_workspace
    I_dcnorm_step = t_sam / t_dc * I_dc
    print(I_dcnorm_step)
    assert False
    #assert np.allclose(ws.extractY().ravel(), I_dcnorm_step.ravel())


if __name__ == '__main__':
    pytest.main()
