# local imports
from drtsans.beam_finder import _calculate_neutron_drop

# third party imports
from numpy.testing import assert_almost_equal
import pytest

# standard imports


def test_calculate_neutron_drop():
    path_length = 15.0  # meters
    wavelength = 18.0  # Angstrom
    gravity_drop = _calculate_neutron_drop(path_length, wavelength)
    assert_almost_equal(gravity_drop, 0.02284, decimal=4)


if __name__ == "__main__":
    pytest.main([__file__])
