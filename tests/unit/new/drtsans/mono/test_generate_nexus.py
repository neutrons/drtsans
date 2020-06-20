import pytest
from drtsans.generate_event_nexus import EventNeXusWriter
from drtsans.generate_event_nexus import convert_histogram_to_events


def test_imports():
    assert EventNeXusWriter


def test_convert_histogram_to_events():
    """

    Returns
    -------

    """
    assert convert_histogram_to_events


if __name__ == '__main__':
    pytest.main(__file__)
