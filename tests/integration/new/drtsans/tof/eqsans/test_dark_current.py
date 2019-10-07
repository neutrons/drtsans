import pytest
import numpy as np
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans.dark_current import normalise_to_workspace

# Data for test work space was from https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/174
x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
# total Dark current intensity (summed over wavelength)
I_dc_total = (x + 5. * y) * 10. + 100
# time duration of a frame
t_frame = 1./60.*1000000.  # in microseconds
# sample measurement tof range
t_low = 500  # us
t_high = 2000  # us
# measurement duration
t_sam = 100  # sample measure for 100s
t_dc = 3600  # dark current collected for 1 hr.
# wavelength range
l_max = 6  # assume wavelength max is 6
l_min = 2.5  # assume wavelength min is 2.5
l_step = 0.1  # assume wavelength binning size is 0.1
# DC per wavelength data
nqbins = int((l_max-l_min+l_step/2.)//l_step)
I_dc_oneimage = I_dc_total / nqbins
# DC data 3D (image + wavelength)
I_dc = np.repeat(I_dc_oneimage, nqbins)
I_dc.shape = 5, 5, -1
#
np.set_printoptions(precision=4, suppress=True)


@pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 5, 'Ny': 5}], indirect=True)
def test_normalize_dark_current(workspace_with_instrument):
    """Test of dark current normalization
    dev - Steven Hahn <hahnse@ornl.gov> Jiao Lin <linjiao@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>

    For details see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/156
    and also https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/174
    """
    # create a workspace with fake events with the given wavelength
    wavelength_bin = np.arange(l_min, l_max+l_step/2., l_step)
    assert wavelength_bin.size == nqbins+1
    # The DC normalization method `normalise_to_workspace` requires an input sample data workspace
    data_ws = workspace_with_instrument(axis_values=wavelength_bin, intensities=np.zeros((5, 5, nqbins)), view='pixel')
    # The dark current data was already created above
    dark_ws = workspace_with_instrument(axis_values=wavelength_bin, intensities=I_dc, view='pixel')
    # Prepare the fake sample logs so the normalization algorith can work
    data_sample_logs = SampleLogs(data_ws)
    data_sample_logs.insert('tof_frame_width_clipped', t_frame - t_low - t_high)
    data_sample_logs.insert('tof_frame_width', t_frame)
    data_sample_logs.insert('duration', t_sam)
    data_sample_logs.insert('wavelength_min', l_min, unit='Angstrom')
    data_sample_logs.insert('wavelength_max', l_max, unit='Angstrom')
    data_sample_logs.insert('wavelength_lead_min', l_min, unit='Angstrom')
    data_sample_logs.insert('wavelength_lead_max', l_max, unit='Angstrom')
    data_sample_logs.insert('is_frame_skipping', False)
    dark_sample_log = SampleLogs(dark_ws)
    dark_sample_log.insert('tof_frame_width', t_frame)
    dark_sample_log.insert('tof_frame_width_clipped', t_frame - t_low - t_high)
    dark_sample_log.insert('duration', t_dc)
    dark_sample_log.insert('wavelength_min', l_min, unit='Angstrom')
    dark_sample_log.insert('wavelength_max', l_max, unit='Angstrom')
    dark_sample_log.insert('wavelength_lead_min', l_min, unit='Angstrom')
    dark_sample_log.insert('wavelength_lead_max', l_max, unit='Angstrom')
    dark_sample_log.insert('is_frame_skipping', False)
    # This is the result documented in ticket 174 by Changwoo
    expected_result = np.array([
        [0.027, 0.0304,  0.0337, 0.0371, 0.0405],
        [0.0438, 0.0472, 0.0506, 0.054,  0.0573],
        [0.0607, 0.0641, 0.0675, 0.0708, 0.0742],
        [0.0776, 0.081,  0.0843, 0.0877, 0.0911],
        [0.0944, 0.0978, 0.1012, 0.1046, 0.1079],
    ])
    # This implements the formula in ticket 174
    I_dcnorm_step = t_sam / t_dc * ((t_frame - t_low - t_high) / t_frame) * (l_step / (l_max - l_min)) * I_dc_total
    # This is obtained by calling the method id drtsans
    result = normalise_to_workspace(dark_ws, data_ws).extractY()[:, 0].reshape(5, 5) * t_sam
    # All three should match
    np.testing.assert_allclose(result, expected_result, atol=1e-4)
    np.testing.assert_allclose(I_dcnorm_step, expected_result, atol=1e-4)
    return


if __name__ == '__main__':
    pytest.main()
