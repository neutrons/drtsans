import pytest
import numpy as np

x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
I_dc = (x + 5. * y) * 10. + 100
t_frame = 1./60.*1000000.  # in microseconds
t_low = 500  # us
t_high = 2000  # us
t_sam = 100  # sample measure for 100s
t_dc = 3600  # dark current collected for 1 hr.
l_max = 6  # assume wavelength max is 6
l_min = 2.5  # assume wavelength min is 2.5
l_step = 0.1  # assume wavelength binning size is 0.1

#np.set_printoptions(precision=4, suppress=True)
@pytest.mark.parametrize('generic_workspace', [{'Nx': 5, 'Ny': 5}], indirect=True)
def test_dark_current_subtract_dark_current(generic_workspace):
    ws = generic_workspace
    # frame_width_clipped / (frame_width * n_bins * duration) * I_dc(x, y)
    frame_width = (t_frame - t_low - t_high) / t_frame
    #print(frame_width)
    # note: the l_max - l_min should be same size(but different unit) as the t_frame - t_low - t_high
    I_dcnorm_step = t_sam / t_dc * ((t_frame - t_low - t_high) / t_frame) * (l_step / (l_max - l_min)) * I_dc
    print(I_dcnorm_step)
    assert False
    #assert np.allclose(ws.extractY().ravel(), I_dcnorm_step.ravel())


if __name__ == '__main__':
    pytest.main()
