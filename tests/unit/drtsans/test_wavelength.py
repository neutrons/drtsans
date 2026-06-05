import pytest

from drtsans.wavelength import tof, from_tof, BROGLIE_NEUTRON, Wband, Wbands


def test_tof():
    wavelength, distance = 4.0, 10.0
    # Without emission delay, tof = distance / velocity = wavelength * distance / BROGLIE_NEUTRON
    assert tof(wavelength, distance) == pytest.approx(wavelength * distance / BROGLIE_NEUTRON)
    # With an emission_delay callable, the delay is added to the flight time
    t0 = 100.0
    assert tof(wavelength, distance, emission_delay=lambda _: t0) == pytest.approx(
        wavelength * distance / BROGLIE_NEUTRON + t0
    )


def test_from_tof():
    wavelength, distance = 4.0, 10.0
    # Without emission delay, from_tof is the exact inverse of tof
    assert from_tof(tof(wavelength, distance), distance=distance) == pytest.approx(wavelength)
    # With emission_delay, from_tof iteratively inverts tof to recover the original wavelength
    from drtsans.tof.eqsans.correct_frame import emission_delay

    assert from_tof(
        tof(wavelength, distance, emission_delay=emission_delay), distance=distance, emission_delay=emission_delay
    ) == pytest.approx(wavelength, abs=0.001)


class TestWband:
    def test_init(self):
        with pytest.raises(ValueError):
            Wband(-1, 0)
            assert False, 'Should have failed "Wband(-1, 0)"'

    def test_width(self):
        assert Wband(1, 2).width == 1
        assert Wband(1, 1).width == 0

    def test_intersect(self):
        b = Wband(1, 2)
        assert b * Wband(0, 0.5) is None
        assert b * Wband(0, 1) is None
        assert b * Wband(0, 1.5) == Wband(1, 1.5)
        assert b * Wband(0, 2) == b
        assert b * Wband(2, 3) is None
        assert b * Wband(2.5, 3) is None

    def test_almost_equal_same_bands(self):
        b1 = Wband(1, 2)
        b2 = Wband(1, 2)
        assert b1.almost_equal(b2) is True

    def test_almost_equal_different_bands(self):
        b1 = Wband(1, 2)
        b2 = Wband(1.1, 2)
        assert b1.almost_equal(b2, atol=1e-6) is False

    def test_almost_equal_same_min(self):
        b1 = Wband(1, 2)
        b2 = Wband(1, 2.1)
        assert b1.almost_equal(b2, atol=1e-6) is False

    def test_almost_equal_same_max(self):
        b1 = Wband(1, 2)
        b2 = Wband(1.1, 2)
        assert b1.almost_equal(b2, atol=1e-6) is False


class TestWbands:
    def test_init(self):
        # initialize with a Wband object
        ws = Wbands(Wband(1, 2))
        assert len(ws) == 1
        # initialize with multiple arguments
        ws = Wbands(Wband(1, 2), Wband(0, 0.5))
        assert len(ws) == 2
        ref = Wbands(Wband(0, 0.5), Wband(1, 2))
        assert ws == ref
        # initialize from iterable
        ws = Wbands([Wband(1, 2), Wband(0, 0.5)])
        assert ws == ref
        # initialize from Wbands object
        bs = Wbands(ws)
        assert bs == ref
        # Mix object types in initializer
        bs = Wbands(Wband(3, 5), ws)
        assert bs == Wbands(Wband(0, 0.5), Wband(1, 2), Wband(3, 5))

    def test_mul(self):
        ws = Wbands(Wband(1, 2), Wband(3, 5))
        # Product of one Wband with one Wbands
        assert ws * Wband(0, 4) == Wbands(Wband(1, 2), Wband(3, 4))
        assert Wband(0, 4) * ws == Wbands(Wband(1, 2), Wband(3, 4))
        assert ws * Wband(1.5, 3.5) == Wbands(Wband(1.5, 2), Wband(3, 3.5))
        assert Wband(1.5, 3.5) * ws == Wbands(Wband(1.5, 2), Wband(3, 3.5))
        # Product of two Wbands
        vs = Wbands(Wband(0, 1.5), Wband(1.7, 4))
        assert ws * vs == Wbands(Wband(1.7, 2), Wband(3, 4), Wband(1, 1.5))
        assert vs * ws == Wbands(Wband(1, 1.5), Wband(3, 4), Wband(1.7, 2))
        # Product of two Wband and one Wbands
        assert Wband(2, 3.5) * ws * Wband(0, 4) == Wbands(Wband(3, 3.5))
        # Product of three Wbands
        intersection = Wbands(Wband(0, 1.5), Wband(2, 3.5)) * ws * vs
        assert intersection == Wbands(Wband(1, 1.5), Wband(3, 3.5))

    def test_getitem(self):
        ws = Wbands(Wband(1, 1.5), Wband(3, 4), Wband(1.7, 2))
        assert ws[1] == Wband(1.7, 2)

    def test_almost_equal_same_bands(self):
        b1 = Wbands(Wband(1, 2), Wband(3, 4))
        b2 = Wbands(Wband(1, 2), Wband(3, 4))
        assert b1.almost_equal(b2) is True

    def test_almost_equal_different_bands(self):
        b1 = Wbands(Wband(1, 2), Wband(3, 4))
        b2 = Wbands(Wband(1, 2.1), Wband(3, 4))
        assert b1.almost_equal(b2, atol=1e-6) is False

    def test_almost_equal_different_length(self):
        b1 = Wbands(Wband(1, 2), Wband(3, 4))
        b2 = Wbands(Wband(1, 2), Wband(3, 4), Wband(5, 6))
        assert b1.almost_equal(b2) is False


if __name__ == "__main__":
    pytest.main([__file__])
