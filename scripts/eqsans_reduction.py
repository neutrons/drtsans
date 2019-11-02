import json
import os
import sys
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
from drtsans.tof import eqsans  # noqa E402
from drtsans.tof.eqsans.convert_to_q import convert_to_q, split_by_frame  # noqa E402
from drtsans.iq import bin_iq_into_linear_q1d, BinningMethod  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D  # noqa E402


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
    config["mask"] = json_conf["maskFileName"]
    config["flux"] = json_conf["beamFluxFileName"]
    config["sensitivity_file_path"] = json_conf["sensitivityFileName"]
    config["bin_width"] = json_conf["wavelenStep"]
    config["low_tof_clip"] = float(json_conf["TOFmin"])
    config["high_tof_clip"] = float(json_conf["TOFmax"])
    config["detector_offset"] = float(json_conf["detectorOffset"])
    config["sample_offset"] = float(json_conf["sampleOffset"])
    config["flux_method"] = "proton charge"  # TODO: what's this?

    # find the beam center
    empty_run = json_params["empty"]["runNumber"]
    empty_fn = "EQSANS_{}".format(empty_run)
    # TODO apply empty flag?
    if empty_run != "":
        db_ws = eqsans.load_events(empty_fn)
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
                                              radius_unit="mm")
        ws = eqsans.apply_transmission_correction(ws,
                                                  trans_workspace=tr_ws)

    # background
    bkg_run = json_params["background"]["runNumber"]
    if bkg_run != "":
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
                                                  radius_unit="mm")
            ws_bck = eqsans.apply_transmission_correction(ws_bck,
                                                          trans_workspace=tr_ws)

        ws = eqsans.subtract_background(ws, background=ws_bck)

    ws /= sample_thickness
    ws *= absolute_scale

    # convert to momentum transfer and split by frame
    all_I_of_q = convert_to_q(ws, mode="scalar")
    I_of_q_by_frame = split_by_frame(ws, all_I_of_q)

    # do the binning and save the files
    numQBins = int(json_conf["numQBins"])
    linear_binning = json_conf["QbinType"] == "linear"
    fig, ax = plt.subplots()
    label = None
    save_suffix = ''
    title = f"reduction log {output_file}"
    for frame_number, result in enumerate(I_of_q_by_frame):
        if linear_binning:
            binned_i_of_q = bin_iq_into_linear_q1d(result.intensity,
                                                   result.error,
                                                   result.mod_q,
                                                   result.delta_q,
                                                   bins=numQBins,
                                                   bin_method=BinningMethod.WEIGHTED)
            if len(I_of_q_by_frame) > 1:
                label = f"frame_{frame_number+1}"
                save_suffix = f"_frame_{frame_number+1}"
            ax.loglog(binned_i_of_q.q, binned_i_of_q.i, label=label)
            filename = os.path.join(json_conf["outputDir"],
                                    json_params["outputFilename"] + save_suffix + '_Iq.txt')
            save_ascii_binned_1D(filename, title, binned_i_of_q)

    ax.legend()
    ax.set_ylabel("Intensity")
    ax.set_xlabel("$Q (\AA^{-1})$")  # noqa W605

    suffix = ".png"
    picture_file = os.path.join(json_conf["outputDir"],
                                json_params["outputFilename"] + suffix)
    fig.savefig(picture_file)
