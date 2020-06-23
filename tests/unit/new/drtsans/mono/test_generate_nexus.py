import pytest
import os
import h5py
from drtsans.h5_buffer import parse_h5_entry
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


def test_create_instrument_node(reference_dir):
    """Test to create an instrument node

    Returns
    -------

    """
    # Parse NeXus file manually
    source_nexus = os.path.join(reference_dir.new.gpsans, 'CG2_9166.nxs.h5')

    # Parse the source HDF5
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = parse_h5_entry(nexus_h5)

    # IDF in XML
    xml_idf = nexus_h5['entry']['instrument']['instrument_xml']['data'][0]

    # Create new instrument node
    test_node = InstrumentNode()
    test_node.set_idf(xml_idf, idf_type=b'text/xml', description=b'XML contents of the instrument IDF')
    test_node.set_instrument_info(target_station_number=1, beam_line=b'CG2', name=b'CG2', short_name=b'CG2')

    # Verify
    # attributes
    source_instrument = source_root.get_child('/entry').get_child('/entry/instrument')

    # attributes
    # cannot get b'NXinstrument'
    assert source_instrument.attributes == test_node.attributes, '{} shall be same as {}' \
                                                                 ''.format(source_instrument.attributes,
                                                                           test_node.attributes)

    # beam line
    for child_name in ['beamline', 'instrument_xml', 'name', 'target_station_number']:
        child_name = f'/entry/instrument/{child_name}'
        src_child_node = source_instrument.get_child(child_name)
        test_child_node = test_node.get_child(child_name)
        src_child_node.match(test_child_node)


if __name__ == '__main__':
    pytest.main(__file__)
