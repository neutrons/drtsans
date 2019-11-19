import json
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
from drtsans.tof import eqsans  # noqa E402
from drtsans.iq import bin_intensity_into_q1d, BinningMethod, bin_iq_into_linear_q2d, BinningParams  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402

debug_mode = False
# TODO: move it to read from a configuration file
default_mask = [{'Bank': '2, 26'}, {'Bank': '19', 'Tube': '2'}, {'Pixel': '1-8, 249-256'}]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("reduction code requires a parameter json string")

    json_string = " ".join(sys.argv[1:])
    json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    output_file = json_params["outputFilename"]

    # set up the configuration
    config = dict()
    json_conf = json_params["configuration"]
    config["mask"] = None
    w = msapi.LoadEmptyInstrument(InstrumentName='EQ-SANS', OutputWorkspace=uwd())
    if json_conf["useDefaultMask"]:
        for d in default_mask:
            msapi.MaskBTP(Workspace=w, **d)
        config["mask"] = msapi.ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace

    config["flux"] = json_conf["beamFluxFileName"]
    config["sensitivity_file_path"] = json_conf["sensitivityFileName"]
    config["dark_current"] = json_conf["darkFileName"]
    config["bin_width"] = json_conf["wavelenStep"]
    default_tof_cut_low = 500
    default_tof_cut_high = 2000
    if json_conf["useTOFcuts"]:
        config["low_tof_clip"] = float(json_conf["TOFmin"])
        config["high_tof_clip"] = float(json_conf["TOFmax"])
    else:
        config["low_tof_clip"] = default_tof_cut_low
        config["high_tof_clip"] = default_tof_cut_high
    config["detector_offset"] = float(json_conf["detectorOffset"])
    config["sample_offset"] = float(json_conf["sampleOffset"])
    config["flux_method"] = "proton charge"  # TODO: what's this?
    if json_conf["useMaskBackTubes"]:
        config["mask_panel"] = "back"

    # find the beam center
    empty_run = json_params["empty"]["runNumber"]
    empty_fn = "EQSANS_{}".format(empty_run)
    # TODO apply empty flag?
    if empty_run != "":
        db_ws = eqsans.load_events(empty_fn)
        if json_conf["useDefaultMask"]:
            for d in default_mask:
                msapi.MaskBTP(Workspace=db_ws, **d)
        center = eqsans.find_beam_center(db_ws)
        config["x_center"] = center[0]
        config["y_center"] = center[1]
        msapi.logger.notice("calculated center {}".format(center))
    else:
        config["x_center"] = 0.025239
        config["y_center"] = 0.0170801

    # load and prepare scattering data
    sample_run = json_params["runNumber"]
    sample_file = "EQSANS_{}".format(sample_run)
    if not output_file:
        output_file = sample_file + "_log.hdf5"
    ws = eqsans.prepare_data(sample_file, **config)
    msapi.logger.warning(str(config))
    # TODO check the next two values if empty
    absolute_scale = float(json_conf["absoluteScale"])
    sample_thickness = float(json_params["thickness"])

    # apply transmission
    transmission_run = json_params["transmission"]["runNumber"]
    apply_transmission = transmission_run.strip() is not ''

    if apply_transmission:
        if float(transmission_run) <= 1:
            msapi.logger.notice('...applying transmission correction with fixed value.')
            ws = eqsans.apply_transmission_correction(ws,
                                                      trans_value=float(transmission_run))
        else:
            msapi.logger.notice('...applying transmission correction with transmission file.')
            transmission_fn = "EQSANS_{}".format(transmission_run)
            ws_tr_sample = eqsans.prepare_data(transmission_fn, **config)
            ws_tr_direct = eqsans.prepare_data(empty_fn, **config)
            tr_ws = eqsans.calculate_transmission(ws_tr_sample,
                                                  ws_tr_direct,
                                                  radius=None,
                                                  radius_unit="mm")
            tr_fn = os.path.join(json_conf["outputDir"],
                                 json_params["outputFilename"] + '_trans.txt')
            print('...saving transmission file ' + tr_fn)
            msapi.SaveAscii(tr_ws, Filename=tr_fn)
            ws = eqsans.apply_transmission_correction(ws,
                                                      trans_workspace=tr_ws)
    else:
        print('...no transmission correction is applied')

    # background
    bkg_run = json_params["background"]["runNumber"]
    if bkg_run.strip() is not '':
        print('...applying bkg_subtraction.')
        bkg_fn = "EQSANS_{}".format(bkg_run)
        bkg_trans_run = json_params["background"]["transmission"]["runNumber"]
        bkg_trans_fn = "EQSANS_{}".format(bkg_trans_run)

        ws_bkg = eqsans.prepare_data(bkg_fn, **config)

        # apply transmission background
        if bkg_trans_run.strip() is not '':
            if float(bkg_trans_run) <= 1:
                print('...applying bkg_transmission correction with fixed value.')
                ws_bkg = eqsans.apply_transmission_correction(ws_bkg,
                                                              trans_value=float(bkg_trans_run))
            else:
                print('...applying bkg_transmission correction with transmission file.')
                ws_bkg_trans = eqsans.prepare_data(bkg_trans_fn, **config)
                ws_cal_tr_bkg = eqsans.calculate_transmission(ws_bkg_trans,
                                                              ws_tr_direct,
                                                              radius=None,
                                                              radius_unit="mm")
                cal_tr_bkg_fn = os.path.join(json_conf["outputDir"],
                                             json_params["outputFilename"] + '_bkg_' + bkg_trans_run + '_trans.txt')
                print('...saving bkg_transmission file ' + cal_tr_bkg_fn)
                msapi.SaveAscii(ws_cal_tr_bkg, Filename=cal_tr_bkg_fn)
                ws_bkg = eqsans.apply_transmission_correction(ws_bkg,
                                                              trans_workspace=ws_cal_tr_bkg)
        else:
            print('...no transmission correction is applied to background')

        ws = eqsans.subtract_background(ws, background=ws_bkg)
    else:
        print('...no bkg_subtraction.')

    ws /= sample_thickness
    ws *= absolute_scale

    # apply user mask
    eqsans.apply_mask(ws, mask=json_conf["maskFileName"])

    if debug_mode:
        filenamews = os.path.join(json_conf["outputDir"], json_params["outputFilename"] + '.nxs')
        msapi.SaveNexus(ws, filenamews)

    # convert to momentum transfer and split by frame
    all_I_of_q = eqsans.convert_to_q(ws, mode="scalar")
    all_I_of_qxqy = eqsans.convert_to_q(ws, mode="azimuthal")
    I_of_q_by_frame = eqsans.split_by_frame(ws, all_I_of_q)
    I_of_qxqy_by_frame = eqsans.split_by_frame(ws, all_I_of_qxqy)

    # do the binning and save the files
    numQBins1D = int(json_conf["numQBins"])
    numQBins2D = int(json_conf["numQxQyBins"])
    linear_binning = json_conf["QbinType"] == "linear"

    # 1D
    fig, ax = plt.subplots()
    label = None
    save_suffix = ''
    title = f"reduction log {output_file}"

    for frame_number, result in enumerate(I_of_q_by_frame):
        if linear_binning:
            q_min = np.min(result.mod_q)
            q_max = np.max(result.mod_q)
            binning = BinningParams(q_min, q_max, numQBins1D)
            binned_i_of_q = bin_intensity_into_q1d(result, binning, True, BinningMethod.WEIGHTED)
            if len(I_of_q_by_frame) > 1:
                label = f"frame_{frame_number+1}"
                save_suffix = f"_frame_{frame_number+1}"
            ax.loglog(binned_i_of_q.mod_q, binned_i_of_q.intensity, label=label)
            filename = os.path.join(json_conf["outputDir"],
                                    json_params["outputFilename"] + save_suffix + '_Iq.txt')
            # print(binned_i_of_q)
            save_ascii_binned_1D(filename, title, binned_i_of_q)

    if label:
        ax.legend()
    ax.set_ylabel("Intensity")
    ax.set_xlabel("$Q (\AA^{-1})$")  # noqa W605

    suffix = ".png"
    picture_file = os.path.join(json_conf["outputDir"],
                                json_params["outputFilename"] + suffix)
    fig.savefig(picture_file)

    # 2D
    for frame_number, result in enumerate(I_of_qxqy_by_frame):
        qx_min = np.min(result.qx)
        qx_max = np.max(result.qx)
        binning_x = BinningParams(qx_min, qx_max, numQBins2D)
        qy_min = np.min(result.qy)
        qy_max = np.max(result.qy)
        binning_y = BinningParams(qy_min, qy_max, numQBins2D)

        binned_i_of_qxqy = bin_iq_into_linear_q2d(result, binning_x, binning_y, BinningMethod.WEIGHTED)
        fig, ax = plt.subplots()
        pcm = ax.pcolormesh(binned_i_of_qxqy.qx, binned_i_of_qxqy.qy, binned_i_of_qxqy.intensity,
                            norm=colors.LogNorm())
        fig.colorbar(pcm, ax=ax)

        suffix = "_qxqy{}".format(frame_number+1)
        picture_file = os.path.join(json_conf["outputDir"],
                                    json_params["outputFilename"] + suffix + '.png')
        fig.savefig(picture_file)
        filename = os.path.join(json_conf["outputDir"],
                                json_params["outputFilename"] + save_suffix + '_Iqxqy.txt')
        save_ascii_binned_2D(filename, title, binned_i_of_qxqy)
