import json
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
from drtsans.tof import eqsans  # noqa E402


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise RuntimeError("reduction code requires a parameter json string")

    json_string = ' '.join(sys.argv[1:])
    json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    # set up the configuration
    config = dict()
    json_conf = json_params["configuration"]
    config['mask'] = json_conf["maskFileName"]
    config['flux'] = json_conf["beamFluxFileName"]
    config['sensitivity_file_path'] = json_conf["sensitivityFileName"]

    # find the beam center
    empty_run = json_params["empty"]["runNumber"]
    empty_fn = "EQSANS_{}".format(empty_run)
    # TODO apply empty flag?
    if empty_run != '':
        db_ws = eqsans.load_events(empty_fn)
        center = eqsans.find_beam_center(db_ws)
        config['x_center'] = center[0]
        config['y_center'] = center[1]
        msapi.logger.notice('calculated center {}'.format(center))
    else:
        config['x_center'] = -0.025239
        config['y_center'] = -0.0170801

    # load and prepare scattering data
    sample_run = json_params["runNumber"]
    sample_file = "EQSANS_{}".format(sample_run)
    ws = eqsans.prepare_data(sample_file, **config)
    # TODO check the next two values if empty
    absolute_scale = float(json_conf["absoluteScale"])
    sample_thickness = float(json_params["thickness"])

    # apply transmission
    # TODO check as flag
    apply_transmission = False
    if apply_transmission:
        transmission_run = json_params["transmission"]["runNumber"]
        transmission_fn = "EQSANS_{}".format(transmission_run)
        ws_tr_sample = eqsans.prepare_data(transmission_fn, **config)
        ws_tr_direct = eqsans.prepare_data(empty_fn, **config)
        tr_ws = eqsans.calculate_transmission(ws_tr_sample,
                                              ws_tr_direct,
                                              radius=None,
                                              radius_unit='mm')
        ws = eqsans.apply_transmission_correction(ws,
                                                  trans_workspace=tr_ws)

    # background
    bkg_run = json_params["background"]["runNumber"]
    bkg_fn = "EQSANS_{}".format(bkg_run)
    bkg_trans_run = json_params["background"]["transmission"]["runNumber"]
    bkg__trans_fn = "EQSANS_{}".format(bkg_trans_run)

    ws_bck = eqsans.prepare_data(bkg_fn, **config)

    # apply transmission background
    if apply_transmission:
        ws_tr_back = eqsans.prepare_data(bkg__trans_fn, **config)
        # ws_tr_direct = eqsans.prepare_data("EQSANS_88973", **config)
        tr_ws = eqsans.calculate_transmission(ws_tr_back,
                                              ws_tr_direct,
                                              radius=None,
                                              radius_unit='mm')
        ws_bck = eqsans.apply_transmission_correction(ws_bck,
                                                      trans_workspace=tr_ws)

    ws = eqsans.subtract_background(ws, background=ws_bck)

    ws /= sample_thickness
    ws *= absolute_scale

    # If frame_skipping we will have more than one table workspace
    table_ws_list = eqsans.prepare_momentum_transfer(
        ws, wavelength_binning=[0.5])

    for index, table_ws in enumerate(table_ws_list):

        # TODO check the configuration-numQbins and configuration_QbinType
        numQBins = int(json_conf["numQBins"])
        log_binning = json_conf["QbinType"] == "log"

        iq_ws = eqsans.cal_iq(table_ws, bins=numQBins, log_binning=log_binning)

        suffix = ""
        if len(table_ws_list) > 1:
            suffix = "_frame_{}".format(index+1)
        outfile = os.path.join(json_conf["outputDir"],
                               json_params["outputFilename"] + suffix)

        eqsans.ave_ascii_1d(iq_ws, str(json_params["outputFilename"] + suffix),
                            outfile)
