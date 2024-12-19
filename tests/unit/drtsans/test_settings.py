import pytest
from drtsans.settings import namedtuplefy


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

    assert type(y1) is type(y2)
    assert type(z1) is type(z2)
    assert type(y1) is not type(z1)
    assert type(y2) is not type(z2)


def test_offline():
    print("this tests runs when offline")


if __name__ == "__main__":
    pytest.main([__file__])
