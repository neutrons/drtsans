"""
    BIOSANS reduction script
"""
import json
import os
import sys
import numpy as np
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans  # noqa E402
from drtsans.stitch import stitch_profiles  # noqa E402
from drtsans.plots import plot_IQmod  # noqa E402
from drtsans.mono import biosans as sans  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402

from common_utils import get_Iq, get_Iqxqy, setup_configuration  # noqa E402

INSTRUMENT = 'BIOSANS'


def apply_transmission(ws, transmission_run, empty_run, cfg):
    """
        Apply transmission
    """
    # TODO: there must be a better way to indicate that we are supplying a transmission value
    try:
        is_value = float(transmission_run) <= 1
    except ValueError:
        is_value = False

    if is_value:
        msapi.logger.notice('Applying transmission correction with fixed value.')
        ws = sans.apply_transmission_correction(ws,
                                                trans_value=float(transmission_run))
    else:
        msapi.logger.notice('Applying transmission correction with transmission file.')

        # We need to see the beam, which is on the main detector
        _mask_detector = cfg['mask_detector']
        cfg['mask_detector'] = 'wing_detector'
        ws_tr_sample = sans.prepare_data(transmission_run, **cfg)
        ws_tr_direct = sans.prepare_data(empty_run, **cfg)
        cfg['mask_detector'] = _mask_detector

        # TODO: use the number of pixels around the beam spot
        tr_ws = sans.calculate_transmission(ws_tr_sample,
                                            ws_tr_direct,
                                            radius=None,
                                            radius_unit="mm")
        ws = sans.apply_transmission_correction(ws, trans_workspace=tr_ws)
    return ws


def reduction(json_params, config):
    """
        Perform the whole reduction
    """
    # Load and prepare scattering data
    # all the run numbers are associated with instrument name
    ws = sans.prepare_data(json_params["instrumentName"] + json_params["runNumber"], **config)

    # Transmission
    transmission_run = json_params["transmission"]["runNumber"]
    if transmission_run.strip() is not '':
        transmission_fn = json_params["instrumentName"] + json_params["transmission"]["runNumber"]
        empty_run = json_params["instrumentName"] + json_params["empty"]["runNumber"]
        apply_transmission(ws, transmission_fn, empty_run, config)

    # Background
    bkg_run = json_params["instrumentName"] + json_params["background"]["runNumber"]
    if bkg_run != "":
        ws_bck = sans.prepare_data(bkg_run, **config)

        # Background transmission
        transmission_run = json_params["background"]["transmission"]["runNumber"]
        if transmission_run.strip() is not '':
            transmission_fn = json_params["instrumentName"] + json_params["background"]["transmission"]["runNumber"]
            empty_run = json_params["instrumentName"] + json_params["empty"]["runNumber"]
            apply_transmission(ws_bck, transmission_fn, empty_run, config)

        # Subtract background
        ws = drtsans.subtract_background(ws, background=ws_bck)
        msapi.logger.notice("Background subtracted")

    # Final normalization
    absolute_scale = float(json_params["configuration"]["absoluteScale"])
    sample_thickness = float(json_params["thickness"])
    ws /= sample_thickness
    ws *= absolute_scale

    # Convert the Q
    wing_label = '_wing' if config['is_wing'] else ''
    q_data = sans.convert_to_q(ws, mode='scalar')
    iq_output = get_Iq(q_data, json_params["configuration"]["outputDir"],
                       json_params["outputFilename"], label=wing_label,
                       linear_binning=json_params["configuration"]["QbinType"] == "linear",
                       nbins=int(json_params["configuration"]["numQBins"]))

    q_data = sans.convert_to_q(ws, mode='azimuthal')
    get_Iqxqy(q_data, json_params["configuration"]["outputDir"],
              json_params["outputFilename"], label=wing_label,
              nbins=int(json_params["configuration"]["numQxQyBins"]))
    return iq_output


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("reduction code requires a parameter json string")

    if os.path.isfile(sys.argv[1]):
        print(sys.argv[1])
        with open(sys.argv[1], 'r') as fd:
            json_params = json.load(fd)
    else:
        json_string = " ".join(sys.argv[1:])
        json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    # set up the configuration
    config = setup_configuration(json_params, INSTRUMENT)

    # Find the beam center
    # TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["instrumentName"] + json_params["empty"]["runNumber"]
    if empty_run != "":
        # Load and compute beam center position
        db_ws = sans.load_events(empty_run,
                                 overwrite_instrument=True, output_workspace=uwd())
        msapi.MaskDetectors(db_ws, ComponentList='wing_detector')
        center = sans.find_beam_center(db_ws)

        # Store the center position for later use
        config['center_x'] = center[0]
        config['center_y'] = center[1]
        config['center_y_wing'] = center[2]
        msapi.logger.notice("Calculated center {}".format(center))
    else:
        msapi.logger.warning("WE NEED A WAY TO PASS A BEAM CENTER")

    # TODO: We should be able to process both the main detector and the wing
    # at the same time, but we need to be able to process the sensitivity in two steps.
    # This could be hidden in the API and done automatically.
    config['is_wing'] = False
    config['mask_detector'] = 'wing_detector'
    iq_1 = reduction(json_params, config)

    config['is_wing'] = True
    config['mask_detector'] = 'detector1'
    if json_params['configuration']['useSensitivityFileName']:
        filename = json_params['configuration']['sensitivityFileName'].replace('_flood_', '_flood_wing_')
        config['sensitivity_file_path'] = filename

    iq_2 = reduction(json_params, config)

    # Stitch the main detector and the wing
    overlap = 0.2
    q_start = np.max(iq_1.mod_q) - overlap * (np.max(iq_1.mod_q) - np.min(iq_1.mod_q))
    q_end = overlap * (np.max(iq_2.mod_q) - np.min(iq_2.mod_q)) + np.min(iq_2.mod_q)
    merged_profile = stitch_profiles(profiles=[iq_1, iq_2],
                                     overlaps=[q_start, q_end])
    filename = os.path.join(json_params["configuration"]["outputDir"],
                            json_params['outputFilename'] + '_merged_Iq.png')
    plot_IQmod([merged_profile], filename, backend='mpl')
