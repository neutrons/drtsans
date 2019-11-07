"""
    BIOSANS reduction script
"""
import json
import os
import sys
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans
from drtsans.mono import biosans
from drtsans.iq import BinningMethod, BinningParams  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D  # noqa E402


def load_data(filename, is_wing=False, center_x=None, center_y=None, center_y_wing=None, output_workspace=None):
    """
        Load data. This should be part of the drtsans.mono API.
    """
    # If we didn't get a file path, make sure Mantid can find the file
    if not os.path.isfile(filename):
        filename = "CG3_{}".format(filename)

    ws = msapi.LoadEventNexus(Filename=filename, OutputWorkspace=output_workspace, LoadNexusInstrumentXML=False)
    ws = msapi.HFIRSANS2Wavelength(ws, OutputWorkspace=output_workspace)
    ws = msapi.SetUncertainties(ws, "sqrtOrOne", OutputWorkspace=output_workspace)

    # Mask wing
    if is_wing:
        msapi.MaskDetectors(ws, ComponentList='detector1')
    else:
        msapi.MaskDetectors(ws, ComponentList='wing_detector')
        
    if center_x is not None and center_y is not None and center_y_wing is not None:
        biosans.center_detector(ws, center_x=center_x, center_y=center_y, center_y_wing=center_y_wing)

    return ws

def process_data(ws, cfg):
    # Dark current
    if cfg['useDarkFileName']:
        dark_ws = load_data(cfg['darkFileName'],
                            is_wing=cfg['is_wing'],
                            center_x=cfg['center_x'],
                            center_y=cfg['center_y'],
                            center_y_wing=cfg['center_y_wing'],
                            output_workspace='CG3_dark')
        biosans.subtract_dark_current(ws, dark_ws)
        msapi.logger.notice("Dark current subtracted")

    # Normalization
    if cfg['normalization'] == 'Time':
        biosans.normalize_by_time(ws)
    elif cfg['normalization'] == 'Monitor':
        biosans.normalize_by_monitor(ws)

    # Solid angle
    if cfg['useSolidAngleCorrection']:
        biosans.solid_angle_correction(ws)
        msapi.logger.notice("Solid angle correction performed")

    # Sensitivity
    if cfg['is_wing']:
        msapi.logger.warning("WE NEED TO DEAL WITH WING SENSITIVITY")
    else:
        sensitivity_filename = cfg['sensitivityFileName']
        drtsans.apply_sensitivity_correction(ws, sensitivity_filename=sensitivity_filename)

    # Transmission

    return ws

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
    config = dict()
    config = json_params['configuration']

    #TODO: We need a way to tell if we need to reduce the main detector, the wing, or both
    config['is_wing'] = False

    # Find the beam center
    #TODO: We need a way to pass a pre-calculated beam center
    empty_run = json_params["empty"]["runNumber"]
    if empty_run != "":
        # Load and compute beam center position
        db_ws = load_data(empty_run, is_wing=False, output_workspace='beam')
        center = biosans.find_beam_center(db_ws)

        # Store the center position for later use
        config['center_x'] = center[0]
        config['center_y'] = center[1]
        config['center_y_wing'] = center[2]
        msapi.logger.notice("Calculated center {}".format(center))
    else:
        msapi.logger.warning("WE NEED A WAY TO PASS A BEAM CENTER")

    # Load and prepare scattering data
    ws = load_data(json_params["runNumber"],
                   is_wing=config['is_wing'],
                   center_x=config['center_x'],
                   center_y=config['center_y'],
                   center_y_wing=config['center_y_wing'],
                   output_workspace="__scattering")
    ws = process_data(ws, config)

    # background
    bkg_run = json_params["background"]["runNumber"]
    if bkg_run != "":
        ws_bck = load_data(json_params["background"]["runNumber"],
                           is_wing=config['is_wing'],
                           center_x=config['center_x'],
                           center_y=config['center_y'],
                           center_y_wing=config['center_y_wing'],
                           output_workspace="__background")
        ws_bck = process_data(ws_bck, config)
        ws = drtsans.subtract_background(ws, background=ws_bck)
        msapi.logger.notice("Background subtracted")

    # Final normalization
    absolute_scale = float(config["absoluteScale"])
    sample_thickness = float(json_params["thickness"])
    ws /= sample_thickness
    ws *= absolute_scale

    # Convert the Q
    q_data = drtsans.convert_to_q(ws, mode='scalar')
    
    #TODO: Need to get min/max Q from either drtsans or the frontend
    linear_binning = config["QbinType"] == "linear"
    bin_params = BinningParams(min=0.02, max=.20,
                               bins=int(config["numQBins"]))
    iq_output = biosans.bin_intensity_into_q1d(q_data,
                                               bin_params=bin_params,
                                               linear_binning=linear_binning,
                                               bin_method=BinningMethod.NOWEIGHT)

    wing_label = '_wing' if config['is_wing'] else ''
    output_file = os.path.join(config["outputDir"],
                               json_params["outputFilename"] + wing_label + '_Iq.txt')
    save_ascii_binned_1D(output_file, "I(Q)", iq_output)

    fig, ax = plt.subplots()
    ax.errorbar(iq_output.mod_q, iq_output.intensity, yerr=iq_output.error, label="I(Q)")
    plt.xlabel('Q [1/A]')
    plt.ylabel('Intensity')
    output_file = os.path.join(config["outputDir"],
                               json_params["outputFilename"] + wing_label + '_Iq.png')

    fig.savefig(output_file)