import pytest
from drtsans.mono.generate_event_nexus import EventNeXusWriter
from drtsans.mono.generate_event_nexus import convert_histogram_to_events
from drtsans.mono.generate_event_nexus import InstrumentNode


def test_imports():
    assert EventNeXusWriter


def test_convert_histogram_to_events():
    """

    Returns
    -------

    """
    assert convert_histogram_to_events


def test_create_instrument_node():
    """Test to create an instrument node

    Returns
    -------

    """
    test_node = InstrumentNode()
    assert test_node

    # test_node.set_instrument_info()
    # test_node.set_idf()


if __name__ == '__main__':
    pytest.main(__file__)
