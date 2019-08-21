import pytest
from pytest import approx
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  normalise_by_time)


def test_normalise_by_time(refd):
    w = load_events('EQSANS_68168', data_dir=refd.new.eqsans)
    d = SampleLogs(w).duration.value
    w = transform_to_wavelength(w)
    y, e = sum(w.readY(42)), sum(w.readE(42))
    w = normalise_by_time(w)
    assert (sum(w.readY(42)), sum(w.readE(42))) == approx((y / d, e / d))
    assert SampleLogs(w).normalizing_duration.value == 'duration'
    w.delete()


if __name__ == '__main__':
    pytest.main()
