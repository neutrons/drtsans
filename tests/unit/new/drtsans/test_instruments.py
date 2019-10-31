import pytest

from drtsans.instruments import InstrumentEnumName, instrument_enum_name, is_time_of_flight


def test_instrument_name(serve_events_workspace):
    assert instrument_enum_name('EQ-SANS') == InstrumentEnumName.EQSANS
    assert str(instrument_enum_name('EQ-SANS')) == 'EQSANS'
    assert str(instrument_enum_name('CG3')) == 'BIOSANS'
    with pytest.raises(ValueError):
        instrument_enum_name('nonexistantsansinstrument')
        assert False, 'Should have generated an exception'
    input_workspace = serve_events_workspace('EQSANS_92353')
    assert instrument_enum_name(input_workspace) == InstrumentEnumName.EQSANS


def test_is_time_of_flight(serve_events_workspace):
    for query in ('EQSANS', 'EQ-SANS', serve_events_workspace('EQSANS_92353')):
        assert is_time_of_flight(query) is True
    for query in ('GPSANS', 'CG2', 'BIOSANS', 'CG3'):
        assert is_time_of_flight(query) is False
    

if __name__ == '__main__':
    pytest.main([__file__])
