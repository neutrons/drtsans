"""
    GPSANS: preparing sensitivities script
    This implements https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
"""
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans  # noqa E402
from drtsans.mono import gpsans as sans  # noqa E402
from drtsans.iq import BinningMethod, BinningParams  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402


INSTRUMENT = 'GPSANS'


def setup_configuration(json_params):
    """Extract configuration

    Parameters
    ----------
    json_params : str
        Json parameter dictionary

    Returns
    -------

    """
    flood_files = json_params['flood files']
    direct_beam_files = json_params['direct beam files']

    if len(flood_files) != len(direct_beam_files):
        raise RuntimeError('Number of flood files ({}) shall be equal to number of direct beam files ({})'
                           ''.format(len(flood_files), len(direct_beam_files)))

    return flood_files, direct_beam_files


if __name__ == "__main__":
    """Main
    """
    if len(sys.argv) < 2:
        raise RuntimeError("Prepare sensitivities requires a parameter json string")

    # Parse the input JSON string for configuration
    if os.path.isfile(sys.argv[1]):
        # case: Input is a Json file
        print(sys.argv[1])
        with open(sys.argv[1], 'r') as fd:
            json_params = json.load(fd)
    else:
        # case: Input is a Json string
        json_string = " ".join(sys.argv[1:])
        json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    output_file = json_params['outputFilename']

    # Set up the configuration
    flood_runs, beam_center_runs = setup_configuration(json_params)

    # Prepare files for preparing sensitivities

    # Load data

    # Mask top and bottom detectors

    # Find beam center for each flood file

    # Mask beam center of each flood file

    # Set uncertainties to each file

    # Calculate sensitivities for each file


