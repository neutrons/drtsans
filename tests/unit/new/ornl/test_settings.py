from __future__ import (absolute_import, division, print_function)

import pytest
from ornl.settings import namedtuplefy


def test_namedtuplefy():

    @namedtuplefy
    def foo(x):
        return dict(foo=x)

    @namedtuplefy
    def goo(x):
        return dict(goo=x)

    y1 = foo(42)
    z1 = goo(24)
    y2 = foo(41)
    z2 = goo(21)

    assert type(y1) == type(y2)
    assert type(z1) == type(z2)
    assert type(y1) != type(z1)
    assert type(y2) != type(z2)


if __name__ == '__main__':
    pytest.main()
