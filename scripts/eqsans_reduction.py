import json
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
from ornl.sans.sns import eqsans  # noqa E402


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
    empty_run = json_params["emptyBeam"]["runNumber"]
    empty_fn = "EQSANS_{}".format(empty_run)
    # TODO apply empty flag?
    if empty_run != '':
        db_ws = eqsans.load_events(empty_fn)
        center = eqsans.center_detector(db_ws)
        config['x_center'] = center[0]
        config['y_center'] = center[1]
        msapi.logger.notice('calculated center {}'.format(center))
    else:
        config['x_center'] = -0.025239
        config['y_center'] = -0.0170801

    # load and prepare scattering data
    sample_run = json_params["scattering"]["runNumber"]
    sample_file = "EQSANS_{}".format(sample_run)
    ws = eqsans.prepare_data(sample_file, **config)
    # TODO check the next two values if empty
    absolute_scale = float(json_conf["absoluteScale"])
    sample_thickness = float(json_params["scattering"]["thickness"])

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
    bkg_trans_run = json_params["background"]["transmittedRunNumber"]
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

    # Retrieve information from data, this will eventually not be needed
    frame_skipping = False
    scale_match = 50000
    wl_f1_min = ws.getRun()['wavelength_lead_min'].value
    wl_f1_max = ws.getRun()['wavelength_lead_max'].value
    if 'wavelength_skip_min' in ws.getRun():
        frame_skipping = True
        wl_f2_min = ws.getRun()['wavelength_skip_min'].value
        wl_f2_max = ws.getRun()['wavelength_skip_max'].value

    # first frame
    # TODO 0.1 is parameter?
    eqsans.prepare_momentum_transfer(ws, wavelength_binning=[wl_f1_min,
                                                             0.1,
                                                             wl_f1_max])
    # TODO check the configuration-numQbins and configuration_QbinType
    numQBins = int(json_params["binning"]["size"])
    log_binning = json_params["binning"]["useLog"]
    eqsans.iq(ws, bins=numQBins, log_binning=log_binning)
    outfile = os.path.join(json_conf["outputDir"],
                           json_params["outputFilename"]+'_frame1.txt')
    msapi.SaveAscii(ws.name() + "_iq",
                    Filename=outfile,
                    WriteSpectrumID=False,
                    WriteXError=True)
    # second frame
    if frame_skipping:
        eqsans.prepare_momentum_transfer(ws,
                                         wavelength_binning=[wl_f2_min,
                                                             0.1,
                                                             wl_f2_max])
        # TODO can we have different binning
        eqsans.iq(ws, bins=numQBins, log_binning=log_binning)
        outfile = os.path.join(json_conf["outputDir"],
                               json_params["outputFilename"]+'_frame2.txt')
        msapi.SaveAscii(ws.name() + "_iq",
                        Filename=outfile,
                        WriteSpectrumID=False,
                        WriteXError=True)
