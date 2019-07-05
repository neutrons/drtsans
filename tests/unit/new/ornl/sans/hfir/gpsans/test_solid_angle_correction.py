import pytest
import ornl.sans.hfir.gpsans as gpsans

# most tests in tests/unit/new/ornl/sans/test_solid_angle.py
# just verify that the function is present here.


def test_load_events():
    assert "solid_angle_correction" in dir(gpsans)


if __name__ == '__main__':
    pytest.main()
