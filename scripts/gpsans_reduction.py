"""
    GPSANS reduction script
"""
from datetime import datetime
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

        transmission_dict = {'value': float(transmission_run),
                             'error': ''}

    else:
        msapi.logger.notice('Applying transmission correction with transmission file.')

        # We need to see the beam, which is on the main detector
        ws_tr_sample = sans.prepare_data(transmission_run, output_suffix='_trans_sample', **transmission_cfg)
        ws_tr_direct = sans.prepare_data(empty_run, output_suffix='_trans_direct', **transmission_cfg)

        tr_ws = sans.calculate_transmission(ws_tr_sample,
                                            ws_tr_direct,
                                            radius=cfg['transmission_radius'],
                                            radius_unit="mm")

        transmission_dict = {'value': tr_ws.extractY(),
                             'error': tr_ws.extractE()}

        # remove the temporary workspaces
        ws_tr_sample.delete()
        ws_tr_direct.delete()

        # apply the calculation
        ws = sans.apply_transmission_correction(ws, trans_workspace=tr_ws)

        # remove transmission correction
        if str(tr_ws) in msapi.mtd:  # protect against non-workspaces
            tr_ws.delete()

    return ws, transmission_dict


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
    file_name = json_params["runNumber"]
    if not os.path.exists(file_name):
        file_name = json_params["instrumentName"] + "_" + file_name
    ws = sans.prepare_data(file_name, output_suffix='_data', **config)

    # Transmission
    transmission_run = json_params["transmission"]["runNumber"]
    sample_transmission_dict = {}
    if transmission_run.strip() != '':
        if not os.path.exists(transmission_run):
            transmission_run = json_params["instrumentName"] + "_" + transmission_run
        empty_run_fn = json_params["empty"]["runNumber"]
        if not os.path.exists(empty_run_fn):
            empty_run_fn = json_params["instrumentName"] + "_" + empty_run_fn
        ws, sample_transmission_dict = apply_transmission(ws, transmission_run, empty_run_fn, config)

    # Background
    bkg_run = json_params["background"]["runNumber"]
    background_transmission_dict = {}
    if bkg_run != '':
        if not os.path.exists(bkg_run):
            bkg_run = json_params["instrumentName"] + "_" + bkg_run
        ws_bck = sans.prepare_data(bkg_run, output_suffix='_bkg', **config)

        # Background transmission
        transmission_run = json_params["background"]["transmission"]["runNumber"]
        if transmission_run.strip() != '':
            if not os.path.exists(transmission_run):
                transmission_run = json_params["instrumentName"] + "_" + transmission_run
            empty_run = json_params["empty"]["runNumber"]
            if not os.path.exists(empty_run):
                empty_run = json_params["instrumentName"] + "_" + empty_run
            ws_bck, background_transmission_dict = apply_transmission(ws_bck, transmission_fn, empty_run, config)

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
    Iq = get_Iq(q_data, json_params["configuration"]["outputDir"],
                json_params["outputFilename"],
                linear_binning=json_params["configuration"]["QbinType"] == "linear",
                weighting=flag_weighted,
                nbins=int(json_params["configuration"]["numQBins"]))

    q_data = sans.convert_to_q(ws, mode='azimuthal')
    Iqxqy = get_Iqxqy(q_data, json_params["configuration"]["outputDir"],
                      json_params["outputFilename"],
                      weighting=flag_weighted,
                      nbins=int(json_params["configuration"]["numQxQyBins"]))

    return {'iq': Iq,
            'iqxqy': Iqxqy,
            'sample_transmission': sample_transmission_dict,
            'background_transmission': background_transmission_dict}


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
    log_json_params = copy.deepcopy(json_params)

    msapi.logger.notice(json.dumps(json_params, indent=2))
    msapi.logger.notice("drtsans version: {}".format(drtsans.__version__))

    output_file = json_params['outputFilename']

    # set up the configuration
    config = setup_configuration(json_params, INSTRUMENT)

    # Find the beam center
    # TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["empty"]["runNumber"]
    if not os.path.exists(empty_run):
        transmission_fn = json_params["instrumentName"] + "_" + empty_run
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

    reduction_dict = reduction(json_params, config)
    Iq = reduction_dict['iq']
    Iqxqy = reduction_dict['iqxqy']
    sample_transmission_dict = reduction_dict['sample_transmission']
    background_transmission_dict = reduction_dict['background_transmission']

    # list of arguments for log file ========================================================
    filename = os.path.join(json_params["configuration"]["outputDir"], output_file + '_reduction_log.hdf')
    starttime = datetime.now().isoformat()
    pythonfile = __file__
    reductionparams = log_json_params
    specialparameters = {'beam_center': {'x': config['center_x'],
                                         'y': config['center_y']},
                         'sample_transmission': sample_transmission_dict,
                         'background_transmission': background_transmission_dict,
                         }
    detectordata = {'main': {'iq': Iq, 'iqxqy': Iqxqy}}
    drtsans.savereductionlog(filename=filename,
                             detectordata=detectordata,
                             reductionparams=reductionparams,
                             pythonfile=pythonfile,
                             starttime=starttime,
                             specialparameters=specialparameters,
                             )
