import pytest

from drtsans import geometry as geo


def test_instrument_name(serve_events_workspace):
    assert geo.instrument_enum_name('EQ-SANS') == geo.InstrumentEnumName.EQSANS
    assert str(geo.instrument_enum_name('EQ-SANS')) == 'EQSANS'
    assert str(geo.instrument_enum_name('CG3')) == 'BIOSANS'
    with pytest.raises(ValueError):
        geo.instrument_enum_name('nonexistantsansinstrument')
        assert False, 'Should have generated an exception'
    input_workspace = serve_events_workspace('EQSANS_92353')
    assert geo.instrument_enum_name(input_workspace) == geo.InstrumentEnumName.EQSANS


if __name__ == '__main__':
    pytest.main([__file__])
