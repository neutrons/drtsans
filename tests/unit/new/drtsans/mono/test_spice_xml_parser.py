import pytest
import os
from drtsans.mono.spice_xml_parser import SpiceXMLParser


def test_get_das_logs(reference_dir):
    """

    Returns
    -------

    """
    # Init the parser
    # test_xml = '/Users/wzz/Projects/SANS/sans-backend/temp/FromXML/CG2_exp315_scan0005_0060.xml'
    test_xml = os.path.join(reference_dir.new.gpsans, 'CG2_exp315_scan0005_0060.xml')
    assert os.path.exists(test_xml), f'Test XML {test_xml} cannot be found'
    xml_parser = SpiceXMLParser(test_xml)

    nodes = ['sample_to_flange', 'sdd', 'source_aperture_size', 'sample_aperture_size',
             'detector_trans', 'source_distance']
    for node_name in nodes:
        xml_node = xml_parser.get_xml_node(node_name)
        print(f'{xml_node.tag} is located by {node_name}')

    for node_name in nodes:
        value, unit = xml_parser.get_node_value(node_name, float)
        print(f'{node_name}: value = {value}, unit = {unit}')

    nexus_spice_log_map = {
        "CG2:CS:SampleToSi": ("sample_to_flange", 'mm'),
        "sample_detector_distance": ("sdd", 'm'),
        "wavelength": ("lambda", 'angstroms'),
        "wavelength_spread": ("dlambda", None),
        "source_aperture_diameter": ("source_aperture_size", 'mm'),
        "sample_aperture_diameter": ("sample_aperture_size", 'mm'),
        "detector_trans_Readback": ("detector_trans", None),
        "source_distance": ("source_distance", 'm'),
        "beamtrap_diameter": ("beamtrap_diameter", 'mm')}

    for nexus_name, spice_tuple in nexus_spice_log_map.items():
        spice_name, default_unit = spice_tuple
        value, unit = xml_parser.get_node_value(spice_name, float)
        if unit is None:
            unit = default_unit
        print(f'{nexus_name}: value = {value}, unit = {unit}')

    xml_parser.read_attenuator()

    xml_parser.close()


if __name__ == '__main__':
    pytest.main(__file__)
