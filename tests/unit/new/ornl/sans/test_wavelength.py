from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_raises

from ornl.sans.wavelength import Wband


class TestWband(object):

    def test_init(self):
        assert_raises(ValueError, Wband, -1, 0)
        assert_raises(ValueError, Wband, 0, -1)
        assert_raises(ValueError, Wband, 1, 0)

    def test_width(self):
        assert Wband(1, 2).width == 1
        assert Wband(1, 1).width == 0

    def test_intersect(self):
        b = Wband(1, 2)
        assert b.intersect(Wband(0, 0.5)) is None
        assert b.intersect(Wband(0, 1)) is None
        assert b.intersect(Wband(0, 1.5)) == Wband(1, 1.5)
        assert b.intersect(Wband(0, 2)) == b
        assert b.intersect(Wband(2, 3)) is None
        assert b.intersect(Wband(2.5, 3)) is None


if __name__ == '__main__':
    pytest.main()
