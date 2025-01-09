import pytest
from packaging.version import parse as parse_version


def test_version():
    from drtsans import __version__

    drtsans_version = __version__
    assert parse_version(drtsans_version) != "unknown"


if __name__ == "__main__":
    pytest.main([__file__])
