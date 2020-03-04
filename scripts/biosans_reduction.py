"""
    BIOSANS reduction script
"""
from datetime import datetime
import json
import os
import sys
import numpy as np
import copy
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans  # noqa E402
from drtsans.stitch import stitch_profiles  # noqa E402
from drtsans.plots import plot_IQmod  # noqa E402
from drtsans.mono import biosans as sans  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D  # noqa E402
from common_utils import get_Iq, get_Iqxqy, setup_configuration  # noqa E402
from drtsans.path import registered_workspace # noqa E402

INSTRUMENT = 'BIOSANS'


def apply_transmission(ws, transmission_run, empty_run, cfg):
    """
        Apply transmission
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
        _mask_detector = cfg['mask_detector']
        cfg['mask_detector'] = 'wing_detector'
        ws_tr_sample = sans.prepare_data(transmission_run, output_suffix=wksp_suffix('_trans_sample', config), **cfg)
        ws_tr_direct = sans.prepare_data(empty_run, output_suffix=wksp_suffix('_trans_direct', config), **cfg)
        cfg['mask_detector'] = _mask_detector

        tr_ws = sans.calculate_transmission(ws_tr_sample,
                                            ws_tr_direct,
                                            radius=cfg['transmission_radius'],
                                            radius_unit="mm")

        transmission_dict = {'value': tr_ws.extractY(),
                             'error': tr_ws.extractE()}

        ws = sans.apply_transmission_correction(ws, trans_workspace=tr_ws)

        # remove transmission correction
        if str(tr_ws) in msapi.mtd:  # protect against non-workspaces
            tr_ws.delete()

    return ws, transmission_dict


def wksp_suffix(suffix, config):
    if config['is_wing']:
        return suffix + '_wing'
    else:
        return suffix


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
    # all the run numbers are associated with instrument name
    file_name = json_params["runNumber"]
    if not os.path.exists(file_name):
        file_name = json_params["instrumentName"] + "_" + file_name
    ws = sans.prepare_data(file_name, output_suffix=wksp_suffix('_data', config), **config)

    # Transmission
    transmission_run = json_params["transmission"]["runNumber"]
    sample_transmission_dict = {}
    if transmission_run.strip() != '':
        if not os.path.exists(transmission_run):
            transmission_run = json_params["instrumentName"] + "_" + transmission_run
        empty_run = json_params["empty"]["runNumber"]
        if not os.path.exists(empty_run):
            empty_run = json_params["instrumentName"] + "_" + empty_run
        ws, sample_transmission_dict = apply_transmission(ws, transmission_run, empty_run, config)

    # Background
    bkg_run = json_params["background"]["runNumber"]
    background_transmission_dict = {}
    if bkg_run != "":
        if os.path.exists(bkg_run):
            bkg_run = json_params["instrumentName"] + "_" + bkg_run
        ws_bck = sans.prepare_data(bkg_run, output_suffix=wksp_suffix('_bkg', config), **config)

        # Background transmission
        transmission_run = json_params["background"]["transmission"]["runNumber"]
        if transmission_run.strip() != '':
            transmission_fn = transmission_run
            empty_run = json_params["empty"]["runNumber"]
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
    wing_label = '_wing' if config['is_wing'] else ''
    q_data = sans.convert_to_q(ws, mode='scalar')
    iq_output = get_Iq(q_data, json_params["configuration"]["outputDir"],
                       json_params["outputFilename"], label=wing_label,
                       linear_binning=json_params["configuration"]["QbinType"] == "linear",
                       weighting=flag_weighted,
                       nbins=int(json_params["configuration"]["numQBins"]))

    q_data = sans.convert_to_q(ws, mode='azimuthal')
    iqxqy_output = get_Iqxqy(q_data, json_params["configuration"]["outputDir"],
                             json_params["outputFilename"], label=wing_label,
                             weighting=flag_weighted,
                             nbins=int(json_params["configuration"]["numQxQyBins"]))

    return {'iq': iq_output,
            'iqxqy': iqxqy_output,
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
    msapi.logger.notice(json.dumps(json_params, indent=2))
    msapi.logger.notice("drtsans version: {}".format(drtsans.__version__))
    log_json_params = copy.deepcopy(json_params)

    # set up the configuration
    config = setup_configuration(json_params, INSTRUMENT)

    # Find the beam center
    # TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["empty"]["runNumber"]
    if empty_run != "":
        if not os.path.exists(empty_run):
            empty_run = json_params["instrumentName"] + "_" + empty_run
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
    reduction_1_dict = reduction(json_params, config)
    iq_1 = reduction_1_dict['iq']
    iqxqy_1 = reduction_1_dict['iqxqy']
    sample_transmission_1 = reduction_1_dict['sample_transmission']
    background_transmission_1 = reduction_1_dict['background_transmission']

    config['is_wing'] = True
    config['mask_detector'] = 'detector1'
    if json_params['configuration']['useSensitivityFileName']:
        filename = json_params['configuration']['sensitivityFileName'].replace('_flood_', '_flood_wing_')
        config['sensitivity_file_path'] = filename

    reduction_2_dict = reduction(json_params, config)
    iq_2 = reduction_2_dict['iq']
    iqxqy_2 = reduction_2_dict['iqxqy']
    sample_transmission_2 = reduction_2_dict['sample_transmission']
    background_transmission_2 = reduction_2_dict['background_transmission']

    # Stitch the main detector and the wing
    overlap = 0.2
    q_start = np.max(iq_1.mod_q) - overlap * (np.max(iq_1.mod_q) - np.min(iq_1.mod_q))
    q_end = overlap * (np.max(iq_2.mod_q) - np.min(iq_2.mod_q)) + np.min(iq_2.mod_q)
    merged_profile = stitch_profiles(profiles=[iq_1, iq_2],
                                     overlaps=[q_start, q_end])
    filename = os.path.join(json_params["configuration"]["outputDir"],
                            json_params['outputFilename'] + '_merged_Iq.png')
    plot_IQmod([merged_profile], filename, backend='mpl')
    filename = os.path.join(json_params["configuration"]["outputDir"],
                            json_params['outputFilename'] + '_merged_Iq.json')
    plot_IQmod([merged_profile], filename, backend='d3')
    filename = os.path.join(json_params["configuration"]["outputDir"],
                            json_params['outputFilename'] + '_merged_Iq.txt')
    save_ascii_binned_1D(filename, "I(Q)", merged_profile)

    # list of arguments for log file =======================================================
    filename = os.path.join(json_params["configuration"]["outputDir"], json_params['outputFilename'] +
                            '_reduction_log.hdf')
    starttime = datetime.now().isoformat()
    # username = 'Neymar'
    pythonfile = __file__
    reductionparams = log_json_params
    specialparameters = {'beam_center': {'x': config['center_x'],
                                         'y': config['center_y'],
                                         'y_wing': config['center_y_wing'],
                                         },
                         'sample_transmission': {'main': sample_transmission_1,
                                                 'wing': sample_transmission_2},
                         'background_transmission': {'main': background_transmission_1,
                                                     'wing': background_transmission_2},
                         }
    detectordata = {'main': {'iq': iq_1, 'iqxqy': iqxqy_1},
                    'wing': {'iq': iq_2, 'iqxqy': iqxqy_2},
                    'combined': {'iq': merged_profile}}
    drtsans.savereductionlog(filename=filename,
                             detectordata=detectordata,
                             reductionparams=reductionparams,
                             pythonfile=pythonfile,
                             starttime=starttime,
                             specialparameters=specialparameters,
                             )
