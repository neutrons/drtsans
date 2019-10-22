import pytest
from pytest import approx
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (load_events, transform_to_wavelength, normalise_by_time)


def test_normalise_by_time(reference_dir):
    w = load_events('EQSANS_68168', data_dir=reference_dir.new.eqsans)
    d = SampleLogs(w).duration.value
    w = transform_to_wavelength(w)
    y, e = sum(w.readY(42)), sum(w.readE(42))
    w = normalise_by_time(w)
    assert (sum(w.readY(42)), sum(w.readE(42))) == approx((y / d, e / d))
    assert SampleLogs(w).normalizing_duration.value == 'duration'
    w.delete()


if __name__ == '__main__':
    pytest.main([__file__])
