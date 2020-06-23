import pytest
import os
import h5py
from drtsans.h5_buffer import parse_h5_entry
from drtsans.mono.generate_event_nexus import EventNeXusWriter
from drtsans.mono.generate_event_nexus import convert_histogram_to_events
from drtsans.mono.generate_event_nexus import InstrumentNode
from drtsans.mono.generate_event_nexus import DasLogNode


def test_imports():
    assert EventNeXusWriter


def test_convert_histogram_to_events():
    """

    Returns
    -------

    """
    assert convert_histogram_to_events


def test_create_das_log_node(reference_dir):
    """Test to create a DAS log node

    Example: /entry/DASlogs/sample_detector_distance

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Parse NeXus file manually for the values nodes
    source_nexus = os.path.join(reference_dir.new.gpsans, 'CG2_9166.nxs.h5')

    # Parse the source HDF5
    nexus_h5 = h5py.File(source_nexus, 'r')
    # source_root = parse_h5_entry(nexus_h5)

    # Get times and value for /entry/DASlogs/sample_detector_distance
    ssd_entry = nexus_h5['entry']['DASlogs']['sample_detector_distance']
    ssd_times = ssd_entry['time'].value
    ssd_start_time = ssd_entry['time'].attrs['start']
    ssd_value = ssd_entry['value'].value
    ssd_value_unit = ssd_entry['value'].attrs['units']

    # Set up a DAS log node
    ssd_test_node = DasLogNode(log_name='/entry/DASlogs/sample_detector_distance',
                               log_times=ssd_times, log_values=ssd_value,
                               start_time=ssd_start_time, log_unit=ssd_value_unit)

    ssd_test_node.set_device_info(device_id=13, device_name=b'Mot-Galil3',
                                  target=b'/entry/DASlogs/CG2:SampleToDetRBV')

    assert ssd_test_node

    # Close HDF5
    nexus_h5.close()


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

    # Close HDF5
    nexus_h5.close()


if __name__ == '__main__':
    pytest.main(__file__)
