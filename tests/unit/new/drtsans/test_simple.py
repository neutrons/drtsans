import pytest
from packaging.version import parse as parse_version


def test_version():
    from drtsans import __version__ as drtsans_version

    assert parse_version("1.0") < parse_version(drtsans_version)


if __name__ == "__main__":
    pytest.main([__file__])
