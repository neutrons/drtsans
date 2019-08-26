import pytest
from pytest import approx
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (load_events, transform_to_wavelength, normalise_by_time)
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


def test_normalization_by_time():
    x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
    z = (x * 0 + y * 0 + 1)  # constant image
    sigma, mu = 0.3, 0  # gaussian image
    d = np.sqrt(x * x + y * y)
    g = np.exp(-((d-mu)**2/(2 * sigma**2)))
    z = (x + 5 * y) * 10 + 100

    I_sam = z  # choose sample data from g or z
    I_sam_err = np.sqrt(z)

    t_sam = 5  # any value
    t_sam_err = 0.2  # 2 percent error

    I_samnorm = I_sam / t_sam
    I_norm_err = np.sqrt((I_sam_err/t_sam)**2 + (I_sam * t_sam_err / t_sam**2)**2)


def test_normalization_by_monitor_spectrum():
    x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
    z = (x*0 + y*0 + 1)  # constant image

    fm = [5, 5, 4, 4, 3, 3, 3, 3, 2, 2]  # flux to monitor ratio
    phi = [20, 40, 30, 25, 20, 10, 5,  5,  5,  5]  # monitor spectrum
    I_sam = []
    for i in range(10):
        I_sam.append(z)  #

    I_samnorm = []
    for i in range(10):
        I_samnorm.append(I_sam[i] / fm[i] / phi[i])


if __name__ == '__main__':
    pytest.main([__file__])
