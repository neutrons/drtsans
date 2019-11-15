"""
    GPSANS: preparing sensitivities script
    This implements https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
"""
import json
import os
import sys
from drtsans.mono.load import load_events
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
import drtsans.mono.gpsans as gp
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


def load_data(data):
    """Load data
    
    :param data: int, str
        Examples: ``55555`` or ``CG3_55555`` or file path.
    :return: 
    """
    ws = load_events(data, output_workspace=None, data_dir=None, overwrite_instrument=False)

    return ws


def mask_detectors(data_ws):
    """Mask detectors in a workspace

    Mask (1) beam center, (2) top and (3) bottom

    :param data_ws:
    :return:
    """
    # Mask beam centers
    circular_ids = circular_mask_from_beam_center(data_ws, radius=3.1415926, unit='mm')
    data_ws = apply_mask(data_ws, mask=circular_ids, panel=None, output_workspace=None)

    # Mask top

    # Mask bottom

    return


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


def center_detector(ws):
    """Center detector to detector beam center

    Parameters
    ----------
    ws : ~mantid.Workspace

    Returns
    -------
    ~mantid.Workspace

    """
    xc, yc = gp.find_beam_center(ws)

    # center the detector
    gp.center_detector(ws, xc, yc)

    return ws


def main(argv):
    """Main function

    Parameters
    ----------
    argv

    Returns
    -------

    """
    if len(argv) < 2:
        raise RuntimeError("Prepare sensitivities requires a parameter json string")

    # Parse the input JSON string for configuration
    if os.path.isfile(argv[1]):
        # case: Input is a Json file
        print(argv[1])
        with open(argv[1], 'r') as fd:
            json_params = json.load(fd)
    else:
        # case: Input is a Json string
        json_string = " ".join(argv[1:])
        json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    output_file = json_params['output_file_name']

    # Set up the configuration
    flood_runs, beam_center_runs = setup_configuration(json_params)

    # Prepare files for preparing sensitivities

    # Load data
    for data in []:
        ws = load_data(data)

    # Mask top and bottom detectors
    mask_detectors([])

    # Find beam center for each flood file

    # Mask beam center of each flood file

    # Set uncertainties to each file

    # Calculate sensitivities for each file


def generate_test_json():
    """Generate test case

    From the commission runs in Nov. 2019 provided by instrument scientist

    585 Flood 6m center 2019/10/31 20:04:31 EDT 2:00:00 6.08 5.36 8
    586 db for 6m center 2019/10/31 22:06:11 EDT 0:02:00 6.08 5.36 8
    587 Flood 6 m 200 off center 2019/10/31 22:12:28 EDT 2:00:00 6.08 5.36 8
    588 db 6 m 200 off center 2019/11/01 00:14:08 EDT 0:02:00 6.08 5.36 8
    589 Flood 6 m 400 off center 2019/11/01 00:20:27 EDT 2:00:09 6.08 5.36 8
    590 db 6 m 400 off center 2019/11/01 02:22:07 EDT 0:02:00 6.08 5.36 8

    Returns
    -------
    str
        A standard JSON string (dict)

    """
    # Set up the dictionary for JSON
    json_dict = {'flood files': [585, 587, 589],
                 'direct beam files': [586, 588, 590],
                 'output_file_name': 'test_gp_sensitivity.nxs'}
    json_str = json.dumps(json_dict)

    return json_str


if __name__ == "__main__":
    # set up a test case if there is no JSON filled in
    if len(sys.argv) == 1:
        argv = ['', generate_test_json()]
    else:
        argv = sys.argv

    # Call main
    main(argv)
