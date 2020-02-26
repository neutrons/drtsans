"""
    GPSANS reduction script
"""
import json
import os
import sys
import warnings
import copy
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans  # noqa E402
from drtsans.mono import gpsans as sans  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402
from drtsans.path import registered_workspace # noqa #402
from common_utils import get_Iq, get_Iqxqy, setup_configuration  # noqa E402

INSTRUMENT = 'GPSANS'


def apply_transmission(ws, transmission_run, empty_run, cfg):
    """
        Apply transmission

        Note: this doesn't currently work because we are missing DAS logs
    """
    transmission_cfg = copy.deepcopy(cfg)
    for k, v in transmission_cfg.items():
        if 'mask' in k:
            transmission_cfg[k] = None
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
        ws_tr_sample = sans.prepare_data(transmission_run, output_suffix='_trans_sample', **transmission_cfg)
        ws_tr_direct = sans.prepare_data(empty_run, output_suffix='_trans_direct', **transmission_cfg)

        tr_ws = sans.calculate_transmission(ws_tr_sample,
                                            ws_tr_direct,
                                            radius=cfg['transmission_radius'],
                                            radius_unit="mm")
        # remove the temporary workspaces
        ws_tr_sample.delete()
        ws_tr_direct.delete()

        # apply the calculation
        ws = sans.apply_transmission_correction(ws, trans_workspace=tr_ws)

        # remove transmission correction
        if str(tr_ws) in msapi.mtd:  # protect against non-workspaces
            tr_ws.delete()

    return ws


def reduction(json_params, config):
    """
        Perform the whole reduction
    """
    sensitivity_workspace = None
    sensitivity_file_path = config['sensitivity_file_path']
    if sensitivity_file_path is not None:
        config.pop('sensitivity_file_path')
        sensitivity_workspace = uwd()
        drtsans.load_sensitivity_workspace(sensitivity_file_path, sensitivity_workspace)
        config['sensitivity_workspace'] = sensitivity_workspace

    dark_current_workspace = None
    dark_current_file_path = config['dark_current']
    if dark_current_file_path is not None:
        dark_current_workspace = uwd()
        sans.load_dark_current_workspace(dark_current_file_path, dark_current_workspace)
        config['dark_current'] = dark_current_workspace

    # Load and prepare scattering data
    ws = sans.prepare_data(json_params["instrumentName"] + '_' + json_params["runNumber"],
                           output_suffix='_data', **config)

    # Transmission
    transmission_run = json_params["transmission"]["runNumber"]
    if transmission_run.strip() != '':
        transmission_fn = json_params["instrumentName"] + '_' + json_params["transmission"]["runNumber"]
        empty_run = json_params["instrumentName"] + '_' + json_params["empty"]["runNumber"]
        ws = apply_transmission(ws, transmission_fn, empty_run, config)

    # Background
    bkg_run = json_params["background"]["runNumber"]
    if bkg_run != '':
        bkg_run = json_params["instrumentName"] + '_' + bkg_run
        ws_bck = sans.prepare_data(bkg_run, output_suffix='_bkg', **config)

        # Background transmission
        transmission_run = json_params["background"]["transmission"]["runNumber"]
        if transmission_run.strip() != '':
            transmission_fn = json_params["instrumentName"] + "_" + transmission_run
            empty_run = json_params["instrumentName"] + "_" + json_params["empty"]["runNumber"]
            ws_bck = apply_transmission(ws_bck, transmission_fn, empty_run, config)

        # Subtract background
        ws = drtsans.subtract_background(ws, background=ws_bck)
        msapi.logger.notice("Background subtracted")

    if registered_workspace(sensitivity_workspace):
        msapi.DeleteWorkspace(sensitivity_workspace)
    if registered_workspace(dark_current_workspace):
        msapi.DeleteWorkspace(dark_current_workspace)

    # Final normalization
    try:
        absolute_scale = float(json_params["configuration"]["absoluteScale"])
    except ValueError:
        absolute_scale = 1.
    sample_thickness = float(json_params["thickness"])
    ws /= sample_thickness
    ws *= absolute_scale

    # Convert the Q
    flag_weighted = False
    if 'useErrorWeighting' in json_params["configuration"].keys():
        if json_params["configuration"]["useErrorWeighting"] == '':
            flag_weighted = False
        else:
            flag_weighted = json_params["configuration"]["useErrorWeighting"]
    q_data = sans.convert_to_q(ws, mode='scalar')
    get_Iq(q_data, json_params["configuration"]["outputDir"],
           json_params["outputFilename"],
           linear_binning=json_params["configuration"]["QbinType"] == "linear",
           weighting=flag_weighted,
           nbins=int(json_params["configuration"]["numQBins"]))

    q_data = sans.convert_to_q(ws, mode='azimuthal')
    get_Iqxqy(q_data, json_params["configuration"]["outputDir"],
              json_params["outputFilename"],
              weighting=flag_weighted,
              nbins=int(json_params["configuration"]["numQxQyBins"]))


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
    msapi.logger.notice("drtsans version: {}".format(drtsans.__version__))

    output_file = json_params['outputFilename']

    # set up the configuration
    config = setup_configuration(json_params, INSTRUMENT)

    # Find the beam center
    # TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["instrumentName"] + "_" + json_params["empty"]["runNumber"]
    if False and empty_run != "":
        # Load and compute beam center position
        db_ws = sans.load_events(empty_run,
                                 overwrite_instrument=True, output_workspace=uwd())
        center = sans.find_beam_center(db_ws)

        # Store the center position for later use
        config['center_x'] = center[0]
        config['center_y'] = center[1]
        msapi.logger.notice("Calculated center {}".format(center))
    else:
        config['center_x'] = 0
        config['center_y'] = 0

        msapi.logger.warning("WE NEED A WAY TO PASS A BEAM CENTER")

    reduction(json_params, config)
