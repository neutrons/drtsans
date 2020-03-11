"""
Test methods to determine the correct sample and detector positions from meta data and overwriting
"""
import pytest
from mantid.simpleapi import LoadEmptyInstrument
from drtsans.samplelogs import SampleLogs
from drtsans.mono.meta_data import get_sample_detector_offset


@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4,
                           'dx': 0.006, 'dy': 0.004, 'zc': 1.25,
                           'l1': -5.}],
                         indirect=True)
def test_zero_offsets(generic_IDF):
    """Test instrument without offset

    Returns
    -------

    """
    # Generate a generic SANS instrument with detector dimension stated in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/178
    with open(r'/tmp/GenericSANS_Definition.xml', 'w') as tmp:
        tmp.write(generic_IDF)
        tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS',
                             OutputWorkspace='test_integration_roi')

    # Add sample log
    sample_logs = SampleLogs(ws)
    sample_logs.insert('sample-detector-distance', 15.6 * 1E3, unit='mm')
    sample_logs.insert('Generic:CS:SampleToSi', 71, unit='mm')

    # Test method
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'Generic:CS:SampleToSi', 0.071)

    # Verify
    assert sample_offset == 0
    assert detector_offset == 0
