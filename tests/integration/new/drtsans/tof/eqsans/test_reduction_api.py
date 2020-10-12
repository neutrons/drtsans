import pytest
import os
import sys
import json
import tempfile
from drtsans.tof.eqsans import validate_reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402


def test_wavelength_step(reference_dir, cleanfile):

    # Get test JSON file
    test_json_file = os.path.join(reference_dir.new.eqsans, 'whatevername.json')

    # Create output directory
    test_output_dir = tempfile.mkdtemp()
    cleanfile(test_output_dir)

    # Load JSON file to configuration
    with open(sys.argv[1], 'r') as fd:
        input_config = json.load(fd)
        input_config = validate_reduction_parameters(input_config)

    # checking if output directory exists, if it doesn't, creates the folder
    input_config["configuration"]["outputDir"] = test_output_dir
    if not os.path.exists(test_output_dir):
        os.makedirs(test_output_dir)

    loaded = load_all_files(input_config)
    out = reduce_single_configuration(loaded, input_config)

    # verify out
    # blabla...


if __name__ == '__main__':
    pytest.main([__file__])
