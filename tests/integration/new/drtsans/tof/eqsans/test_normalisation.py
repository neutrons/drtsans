import pytest
from pytest import approx
<<<<<<< HEAD:tests/integration/new/drtsans/tof/eqsans/test_normalisation.py
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (load_events, transform_to_wavelength, normalise_by_time)
=======
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  normalise_by_time, normalise_by_monitor)
>>>>>>> eaf0f67... checkpoint:tests/integration/new/ornl/sans/sns/eqsans/test_normalisation.py
import numpy as np


def test_normalise_by_time(reference_dir):
    w = load_events('EQSANS_68168', data_dir=reference_dir.new.eqsans)
    d = SampleLogs(w).duration.value
    w = transform_to_wavelength(w)
    y, e = sum(w.readY(42)), sum(w.readE(42))
    w = normalise_by_time(w)
    assert (sum(w.readY(42)), sum(w.readE(42))) == approx((y / d, e / d))
    assert SampleLogs(w).normalizing_duration.value == 'duration'
    w.delete()


x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
sigma, mu = 0.3, 0  # gaussian image
d = np.sqrt(x * x + y * y)
g = np.exp(-((d - mu) ** 2 / (2 * sigma ** 2)))

@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x,
                           'intensities': g}],
                         indirect=True)
def test_normalization_by_time(generic_workspace):
    ws = generic_workspace
    I_sam = g  # choose sample data from g or z

    t_sam = 5  # any value
    SampleLogs(ws).insert('timer', t_sam, 'Second')
    I_samnorm = I_sam / t_sam
    ws_samnorm = normalise_by_time(ws, log_key='timer')
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


z = (x*0 + y*0 + 1)  # constant image
z_sam = []
for i in range(10):
    z_sam.append(z)


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x,
                           'intensities': z_sam}],
                         indirect=True)
def test_normalization_by_monitor_spectrum(generic_workspace):
    ws = generic_workspace
    fm_ws = ws.clone()
    phi_ws = ws.clone()
    SampleLogs(ws).insert('is_frame_skipping', False)
    fm = [5, 5, 4, 4, 3, 3, 3, 3, 2, 2]  # flux to monitor ratio
    phi = [20, 40, 30, 25, 20, 10, 5,  5,  5,  5]  # monitor spectrum
    I_sam = []
    for i in range(10):
        I_sam.append(z)  #

    I_samnorm = []
    for i in range(10):
        I_samnorm.append(I_sam[i] / fm[i] / phi[i])

    out = normalise_by_monitor(ws,fm_ws,phi_ws)


if __name__ == '__main__':
    pytest.main([__file__])
