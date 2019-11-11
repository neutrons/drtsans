"""
    BIOSANS reduction script
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
from drtsans.mono import biosans as sans  # noqa E402
from drtsans.iq import BinningMethod, BinningParams  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402


INSTRUMENT = 'BIOSANS'


def load_data(filename, is_wing=False, output_workspace=None):
    """
        Load data. This should be part of the drtsans.mono API.
        TODO: Need to deal with detector and sample offsets.
    """
    if output_workspace is None:
        output_workspace = uwd()

    # If we didn't get a file path, make sure Mantid can find the file
    if not os.path.isfile(filename):
        filename = "{}{}".format(INSTRUMENT, filename)

    ws = msapi.LoadEventNexus(Filename=filename, OutputWorkspace=output_workspace, LoadNexusInstrumentXML=False)
    ws = msapi.HFIRSANS2Wavelength(ws, OutputWorkspace=output_workspace)
    ws = msapi.SetUncertainties(ws, "sqrtOrOne", OutputWorkspace=output_workspace)

    if is_wing:
        msapi.MaskDetectors(ws, ComponentList='detector1')
    else:
        msapi.MaskDetectors(ws, ComponentList='wing_detector')

    return ws


def setup_configuration(json_params):
    """
        Extract configuration
    """
    config = json_params['configuration']

    # Masking
    # TODO: get default mask from configuration
    config["mask"] = None
    default_mask = []
    w = msapi.LoadEmptyInstrument(InstrumentName=INSTRUMENT, OutputWorkspace=uwd())
    if config["useDefaultMask"]:
        for d in default_mask:
            msapi.MaskBTP(Workspace=w, **d)
        config["mask"] = msapi.ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace

    if config["useMaskBackTubes"]:
        msapi.logger.notice("Masking back tubes is not implemented")

    return config


def process_data(filename, cfg, output_workspace=None):
    """
        Load and process data.
        This should be in drtsans.mono.biosans
    """
    ws = load_data(filename, is_wing=cfg['is_wing'],
                   output_workspace=output_workspace)

    if cfg['center_x'] is not None and cfg['center_y'] is not None:
        sans.center_detector(ws, center_x=cfg['center_x'], center_y=cfg['center_y'],
                                center_y_wing=cfg['center_y_wing'])

    # Dark current
    if cfg['useDarkFileName']:
        dark_ws = load_data(cfg['darkFileName'],
                            is_wing=cfg['is_wing'],
                            output_workspace=uwd())
        sans.subtract_dark_current(ws, dark_ws)
        msapi.logger.notice("Dark current subtracted")

    # Normalization
    if cfg['normalization'] == 'Time':
        sans.normalize_by_time(ws)
    elif cfg['normalization'] == 'Monitor':
        sans.normalize_by_monitor(ws)

    # Solid angle
    if cfg['useSolidAngleCorrection']:
        sans.solid_angle_correction(ws)
        msapi.logger.notice("Solid angle correction performed")

    # Sensitivity
    sensitivity_filename = cfg['sensitivityFileName']
    if cfg['is_wing']:
        msapi.logger.warning("WE NEED TO DEAL WITH WING SENSITIVITY")
        sensitivity_filename = sensitivity_filename.replace('_flood_', '_flood_wing_')
    if os.path.isfile(sensitivity_filename):
        drtsans.apply_sensitivity_correction(ws, sensitivity_filename=sensitivity_filename)
    else:
        msapi.logger.warning("File doesn't exist: %s" % sensitivity_filename)

    return ws


def apply_transmission(ws, transmission_run, empty_run, cfg):
    """
        Apply transmission

        Note: this doesn't currently work because we are missing DAS logs
    """
    # TODO: there must be a better way to indicate that we are supplying a transmission value
    try:
        is_value = float(transmission_run) <= 1
    except:
        is_value = False

    if is_value:
        msapi.logger.notice('Applying transmission correction with fixed value.')
        ws = sans.apply_transmission_correction(ws,
                                                trans_value=float(transmission_run))
    else:
        msapi.logger.notice('Applying transmission correction with transmission file.')

        # We need to see the beam, which is on the main detector
        _is_wing = cfg['is_wing']
        cfg['is_wing'] = False
        ws_tr_sample = process_data(transmission_run, cfg)
        ws_tr_direct = process_data(empty_run, cfg)
        cfg['is_wing'] = _is_wing

        # TODO: use the number of pixels around the beam spot
        tr_ws = sans.calculate_transmission(ws_tr_sample,
                                            ws_tr_direct,
                                            radius=None,
                                            radius_unit="mm")
        ws = sans.apply_transmission_correction(ws,
                                                trans_workspace=tr_ws)
    return ws


def get_Iq(ws, json_params, config):
    """
        Compute I(q) from corrected workspace
    """
    q_data = drtsans.convert_to_q(ws, mode='scalar')

    # TODO: Need to get min/max Q from either drtsans or the frontend
    linear_binning = config["QbinType"] == "linear"
    q_min = np.min(q_data.mod_q)
    q_max = np.max(q_data.mod_q)
    bin_params = BinningParams(min=q_min, max=q_max,
                               bins=int(config["numQBins"]))
    iq_output = sans.bin_intensity_into_q1d(q_data,
                                            bin_params=bin_params,
                                            linear_binning=linear_binning,
                                            bin_method=BinningMethod.NOWEIGHT)

    wing_label = '_wing' if config['is_wing'] else ''
    output_file = os.path.join(config["outputDir"],
                               json_params["outputFilename"] + wing_label + '_Iq.txt')
    save_ascii_binned_1D(output_file, "I(Q)", iq_output)

    fig, ax = plt.subplots()
    ax.errorbar(iq_output.mod_q, iq_output.intensity, yerr=iq_output.error, label="I(Q)")
    ax.set_xlabel("$Q (\AA^{-1})$")  # noqa W605
    plt.ylabel('Intensity')
    output_file = os.path.join(config["outputDir"],
                               json_params["outputFilename"] + wing_label + '_Iq.png')

    fig.savefig(output_file)


def get_Iqxqy(ws, json_params, config):
    """
        Compute I(qx,qy) from corrected workspace
    """
    q_data = drtsans.convert_to_q(ws, mode='azimuthal')
    qx_min = np.min(q_data.qx)
    qx_max = np.max(q_data.qx)

    binning_x = BinningParams(qx_min, qx_max, int(config["numQxQyBins"]))
    qy_min = np.min(q_data.qy)
    qy_max = np.max(q_data.qy)
    binning_y = BinningParams(qy_min, qy_max, int(config["numQxQyBins"]))

    iq_output = sans.bin_iq_into_linear_q2d(q_data,
                                            qx_bin_params=binning_x,
                                            qy_bin_params=binning_y,
                                            method=BinningMethod.NOWEIGHT)

    wing_label = '_wing' if config['is_wing'] else ''
    filename = os.path.join(config["outputDir"],
                            json_params["outputFilename"] + wing_label + '_Iqxqy.txt')
    save_ascii_binned_2D(filename, "I(Qx,Qy)", iq_output)

    fig, ax = plt.subplots()
    pcm = ax.pcolormesh(iq_output.qx, iq_output.qy, iq_output.intensity,
                        norm=colors.LogNorm())
    fig.colorbar(pcm, ax=ax)
    picture_file = os.path.join(config["outputDir"],
                                json_params["outputFilename"] + wing_label + '_Iqxqy.png')
    fig.savefig(picture_file)


def reduction(json_params, config):
    """
        Perform the whole reduction
    """
    # Load and prepare scattering data
    ws = process_data(json_params["runNumber"], config)

    # Transmission
    transmission_run = json_params["transmission"]["runNumber"]
    if transmission_run.strip() is not '':
        empty_run = json_params["empty"]["runNumber"]
        apply_transmission(ws, transmission_run, empty_run, config)

    # Background
    bkg_run = json_params["background"]["runNumber"]
    if bkg_run != "":
        ws_bck = process_data(json_params["background"]["runNumber"], config)

        # Background transmission
        transmission_run = json_params["background"]["transmission"]["runNumber"]
        if transmission_run.strip() is not '':
            empty_run = json_params["empty"]["runNumber"]
            apply_transmission(ws_bck, transmission_run, empty_run, config)

        # Subtract background
        ws = drtsans.subtract_background(ws, background=ws_bck)
        msapi.logger.notice("Background subtracted")

    # Final normalization
    absolute_scale = float(config["absoluteScale"])
    sample_thickness = float(json_params["thickness"])
    ws /= sample_thickness
    ws *= absolute_scale

    # Apply user mask
    if config['useMaskFileName']:
        drtsans.mask_utils.apply_mask(ws, mask=config["maskFileName"])

    # Convert the Q
    get_Iq(ws, json_params, config)
    get_Iqxqy(ws, json_params, config)


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

    output_file = json_params['outputFilename']

    # set up the configuration
    config = setup_configuration(json_params)

    # Find the beam center
    # TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["empty"]["runNumber"]
    if empty_run != "":
        # Load and compute beam center position
        db_ws = load_data(empty_run, is_wing=False)
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
    reduction(json_params, config)

    config['is_wing'] = True
    reduction(json_params, config)
