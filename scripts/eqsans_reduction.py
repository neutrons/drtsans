from datetime import datetime
import copy
import json
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
import drtsans  # noqa E402
from drtsans.tof import eqsans  # noqa E402
from drtsans.iq import bin_intensity_into_q1d, BinningMethod, bin_intensity_into_q2d  # noqa E402
from drtsans.iq import determine_1d_linear_bins, determine_1d_log_bins  # noqa E402
from drtsans.instruments import extract_run_number # noqa E402
from drtsans.tof.eqsans import cfg  # noqa E402
from drtsans.samplelogs import SampleLogs  # noqa E402
from common_utils import get_Iqxqy  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402
from drtsans.path import registered_workspace # noqa E402
from drtsans.dataobjects import load_iqmod, save_iqmod # noqa E402


def _represents_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def _get_configuration_file_parameters(sample_run):
    try:
        configuration_file_parameters = cfg.load_config(source=sample_run)
    except RuntimeError as e:
        msapi.logger.error(e)
        msapi.logger.warning('Not using previous configuration')
        configuration_file_parameters = {}
    return configuration_file_parameters


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("reduction code requires a parameter json string")
    if os.path.isfile(sys.argv[1]):
        # print(sys.argv[1])
        with open(sys.argv[1], 'r') as fd:
            json_params = json.load(fd)
    else:
        json_string = " ".join(sys.argv[1:])
        json_params = json.loads(json_string)
    log_json_params = copy.deepcopy(json_params)
    msapi.logger.notice(json.dumps(json_params, indent=2))
    msapi.logger.notice("drtsans version: {}".format(drtsans.__version__))

    output_file = json_params["outputFilename"]
    sample_run = json_params["runNumber"]
    if _represents_int(sample_run):
        configuration_file_parameters = _get_configuration_file_parameters(sample_run)
    else:
        configuration_file_parameters = _get_configuration_file_parameters(extract_run_number(sample_run))

    config = dict()
    json_conf = json_params["configuration"]

    config["mask"] = None  # by default, no mask.
    if json_conf["useDefaultMask"]:
        # load masks defined in the eq-sans configuration file
        config["mask"] = configuration_file_parameters['combined mask']
    if json_conf["useMaskBackTubes"]:
        config["mask_panel"] = "back"

    config["detector_offset"] = float(json_conf["detectorOffset"])
    config["sample_offset"] = float(json_conf["sampleOffset"])

    config["bin_width"] = float(json_conf["wavelenStep"])

    # [cd 2/10/2020] weight option added
    flag_weighted = False
    if 'useErrorWeighting' in json_conf.keys():
        if json_conf["useErrorWeighting"] == '':
            flag_weighted = False
        else:
            flag_weighted = json_conf["useErrorWeighting"]
        msapi.logger.notice('...weighting option selected weighted = ' + str(flag_weighted))
    else:
        msapi.logger.notice('...default weighting option used (NOWEIGHT).')

    # [CD 3/2/2020] log binning options
    flag_logqbinsperdecade = None
    if 'LogQBinsPerDecade' in json_conf.keys():
        if json_conf["LogQBinsPerDecade"] != '':
            flag_logqbinsperdecade = json_conf["LogQBinsPerDecade"]
    msapi.logger.notice('...LogQBinsPerDecade : ' + str(flag_logqbinsperdecade))
    flag_logqbinsdecadecenter = False
    if 'LogQBinsDecadeCenter' in json_conf.keys():
        if json_conf["LogQBinsDecadeCenter"] != '':
            flag_logqbinsdecadecenter = json_conf["LogQBinsDecadeCenter"]
    msapi.logger.notice('...LogQBinsDecadeCenter : ' + str(flag_logqbinsdecadecenter))
    flag_logqbinsevendecade = False
    if 'LogQBinsEvenDecade' in json_conf.keys():
        if json_conf["LogQBinsEvenDecade"] != '':
            flag_logqbinsevendecade = json_conf["LogQBinsEvenDecade"]
    msapi.logger.notice('...LogQBinsEvenDecade : ' + str(flag_logqbinsevendecade))

    default_tof_cut_low = 500
    default_tof_cut_high = 2000
    if json_conf["useTOFcuts"]:
        config["low_tof_clip"] = float(json_conf["TOFmin"])
        config["high_tof_clip"] = float(json_conf["TOFmax"])
    else:
        config["low_tof_clip"] = default_tof_cut_low
        config["high_tof_clip"] = default_tof_cut_high

    dark_current_workspace = None
    dark_current_filename = json_conf["darkFileName"]
    if dark_current_filename is not None:
        dark_current_workspace = uwd()
        eqsans.load_dark_current_workspace(dark_current_filename, dark_current_workspace)
        config["dark_current"] = dark_current_workspace

    # [CD, 2/1/2020] flux_method is now taking "normalization" value from the json file.
    # [CD, 2/3/2020] do if statement depending on the normalization method.
    if json_conf["normalization"] == "Monitor":
        config["flux_method"] = 'monitor'
        config["flux"] = json_conf["fluxMonitorRatioFile"]
    elif json_conf["normalization"] == "Time":
        config["flux_method"] = 'time'
        config["flux"] = json_conf["beamFluxFileName"]
    # [CD, 2/3/2020] be careful. Shaman use 'Total charge' but nomalize_by_flux use 'proton charge'
    elif json_conf["normalization"] == "Total charge":
        config["flux_method"] = 'proton charge'
        config["flux"] = json_conf["beamFluxFileName"]
    else:
        config["flux_method"] = None

    if json_conf["useSolidAngleCorrection"] is False:
        config["solid_angle"] = False
    else:
        config["solid_angle"] = True
    logoutput = '...solid angle...' + str(config["solid_angle"])
    msapi.logger.notice(logoutput)

    sensitivity_workspace = None
    sensitivity_file_path = json_conf["sensitivityFileName"]
    if sensitivity_file_path is not None:
        sensitivity_workspace = uwd()
        drtsans.load_sensitivity_workspace(sensitivity_file_path, sensitivity_workspace)
        config['sensitivity_workspace'] = sensitivity_workspace

    # find the beam center
    empty_run = json_params["empty"]["runNumber"]
    empty_fn = json_params["instrumentName"] + '_' + empty_run
    # TODO apply empty flag?
    if empty_run != "":
        db_ws = eqsans.load_events(empty_fn)
        if json_conf["useDefaultMask"]:
            eqsans.apply_mask(db_ws, mask=config['mask'])
        center = eqsans.find_beam_center(db_ws)
        config["center_x"] = center[0]
        config["center_y"] = center[1]
        msapi.logger.notice("calculated center {}".format(center))
    else:
        config["center_x"] = 0.025239
        config["center_y"] = 0.0170801
        msapi.logger.notice("use default center (0.025239, 0.0170801)")
    # load and prepare scattering data
    sample_file = "EQSANS_{}".format(sample_run)
    if not output_file:
        output_file = sample_file

    ws = eqsans.prepare_data(sample_file, output_suffix='_sample', **config)
    msapi.logger.warning(str(config))
    # TODO check the next two values if empty
    # [CD, 1/30/2020] CD has updated this part.
    if json_conf["absoluteScale"] != "":
        absolute_scale = float(json_conf["absoluteScale"])
    else:
        absolute_scale = 1.0

    if json_params["thickness"] != "":
        sample_thickness = float(json_params["thickness"])
    else:
        sample_thickness = 1.0

    # [CD, 1/30/2020] even though the parametername is nPixelsRadius...
    # [CD, 1/30/2020] the value is actually radius in mm unit
    # [CD, 1/30/2020] for future, I suggest changing the parameter name
    # [CD, 1/30/2020] to be "mmRadiusForTransmission" but this needs to be coordinated with UI team
    try:
        rad_trans = float(json_conf["mmRadiusForTransmission"])
    except ValueError:
        rad_trans = None

    # [CD, 1/30/2020] direct beam loading for transmission should happen before apply transmission if statement.
    # [CD, 1/30/2020] otherwise, ws_tr_direct may be missing for bkg transmission calculation
    # [CD, 1/30/2020] current script does not have an option to use different direct beam for the tranmission.
    # [CD, 1/30/2020] empty beam that is used to find beam center is always used as transmission.
    ws_tr_direct = eqsans.prepare_data(empty_fn, output_suffix='_trans_direct', **config)

    # apply transmission
    transmission_run = json_params["transmission"]["runNumber"].strip()
    transmission_value = json_params["transmission"]["value"].strip()
    apply_transmission_run = transmission_run != ''
    if transmission_value != '':
        apply_transmission_value = float(transmission_value) < 1.0
    else:
        apply_transmission_value = False

    sample_transmission_dict = {}
    if apply_transmission_run or apply_transmission_value:
        if apply_transmission_value:
            sample_transmission_dict['value'] = float(transmission_value)
            msapi.logger.notice('...applying transmission correction with fixed value.')
            ws = eqsans.apply_transmission_correction(ws,
                                                      trans_value=float(transmission_value))
        else:
            msapi.logger.notice('...applying transmission correction with transmission file.')
            transmission_fn = "EQSANS_{}".format(transmission_run)
            ws_tr_sample = eqsans.prepare_data(transmission_fn, output_suffix='_trans_sample', **config)
            raw_tr_ws = eqsans.calculate_transmission(ws_tr_sample,
                                                      ws_tr_direct,
                                                      radius=rad_trans,
                                                      radius_unit="mm",
                                                      fit_function='')
            tr_ws = eqsans.calculate_transmission(ws_tr_sample,
                                                  ws_tr_direct,
                                                  radius=rad_trans,
                                                  radius_unit="mm")
            sample_transmission_dict['value'] = tr_ws.extractY()
            sample_transmission_dict['error'] = tr_ws.extractE()

            # [CD, 1/30/2020] radius default input has changed from "None" to "rad_trans"
            # [CD, 1/30/2020] Here, we need both fitted transmission and raw transmission as return values
            # [CD, 1/30/2020] and save both of them as ascii file
            # [CD, 1/30/2020] (or in a single transmission file with multiple columns)
            raw_tr_fn = os.path.join(json_conf["outputDir"],
                                     json_params["outputFilename"] + '_raw_trans.txt')
            msapi.SaveAscii(raw_tr_ws, Filename=raw_tr_fn)
            tr_fn = os.path.join(json_conf["outputDir"],
                                 json_params["outputFilename"] + '_trans.txt')
            msapi.SaveAscii(tr_ws, Filename=tr_fn)
            ws = eqsans.apply_transmission_correction(ws,
                                                      trans_workspace=tr_ws)
    else:
        msapi.logger.warning('...no transmission correction is applied')

    # background
    bkg_run = json_params["background"]["runNumber"]
    background_transmission_dict = {}
    if bkg_run.strip() != '':
        msapi.logger.notice('...applying bkg_subtraction.')
        bkg_fn = "EQSANS_{}".format(bkg_run)

        ws_bkg = eqsans.prepare_data(bkg_fn, output_suffix='_bkg', **config)

        # apply transmission background
        bkg_transmission_run = json_params["background"]["transmission"]["runNumber"].strip()
        bkg_transmission_value = json_params["background"]["transmission"]["value"].strip()
        bkg_apply_transmission_run = bkg_transmission_run != ''
        if bkg_transmission_value != '':
            bkg_apply_transmission_value = float(json_params["transmission"]["value"]) < 1.0
        else:
            bkg_apply_transmission_value = False
        if bkg_apply_transmission_run or bkg_apply_transmission_value:
            if bkg_apply_transmission_value:
                msapi.logger.notice('...applying bkg_transmission correction with fixed value.')
                ws_bkg = eqsans.apply_transmission_correction(ws_bkg,
                                                              trans_value=float(bkg_transmission_value))
                background_transmission_dict['value'] = float(bkg_transmission_value)
            else:
                msapi.logger.notice('...applying bkg_transmission correction with transmission file.')
                bkg_trans_fn = "EQSANS_{}".format(bkg_transmission_run)
                ws_bkg_trans = eqsans.prepare_data(bkg_trans_fn, output_suffix='_bkg_trans', **config)
                ws_cal_raw_tr_bkg = eqsans.calculate_transmission(ws_bkg_trans,
                                                                  ws_tr_direct,
                                                                  radius=rad_trans,
                                                                  radius_unit="mm",
                                                                  fit_function='')
                ws_cal_tr_bkg = eqsans.calculate_transmission(ws_bkg_trans,
                                                              ws_tr_direct,
                                                              radius=rad_trans,
                                                              radius_unit="mm")

                background_transmission_dict['background_raw_transmission'] = {'value': ws_cal_raw_tr_bkg.extractY(),
                                                                               'error': ws_cal_raw_tr_bkg.extractE()}
                background_transmission_dict['background_transmission'] = {'value': ws_cal_tr_bkg.extractY(),
                                                                           'error': ws_cal_tr_bkg.extractE()}

                # [CD, 1/30/2020] radius default input has changed from "None" to "rad_trans"
                # [CD, 1/30/2020] Here, we need both fitted transmission and raw transmission as return values
                # [CD, 1/30/2020] and save both of them as ascii file
                # [CD, 1/30/2020] (or in a single transmission file with multiple columns)
                cal_raw_tr_bkg_fn = os.path.join(json_conf["outputDir"],
                                                 json_params["outputFilename"] +
                                                 '_bkg_' + bkg_transmission_run + '_raw_trans.txt')
                cal_tr_bkg_fn = os.path.join(json_conf["outputDir"],
                                             json_params["outputFilename"] +
                                             '_bkg_' + bkg_transmission_run + '_trans.txt')
                msapi.logger.notice('...saving bkg_transmission file ' + cal_tr_bkg_fn)
                msapi.SaveAscii(ws_cal_tr_bkg, Filename=cal_tr_bkg_fn)
                msapi.SaveAscii(ws_cal_raw_tr_bkg, Filename=cal_raw_tr_bkg_fn)
                ws_bkg = eqsans.apply_transmission_correction(ws_bkg,
                                                              trans_workspace=ws_cal_tr_bkg)
        else:
            msapi.logger.warning('...no transmission correction is applied to background')

        ws = eqsans.subtract_background(ws, background=ws_bkg)
    else:
        msapi.logger.notice('...no bkg_subtraction.')

    if registered_workspace(sensitivity_workspace):
        msapi.DeleteWorkspace(sensitivity_workspace)
    if registered_workspace(dark_current_workspace):
        msapi.DeleteWorkspace(dark_current_workspace)

    ws /= sample_thickness
    ws *= absolute_scale

    # apply user mask
    eqsans.apply_mask(ws, mask=json_conf["maskFileName"])

    # save the nexus file with all corrections
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
    # [AS, 2/4/2010]
    fig, ax = plt.subplots()
    label = None
    save_suffix = ''
    title = f"reduction log {output_file}"

    log_binned_i_of_q = {}
    for frame_number, result in enumerate(I_of_q_by_frame):
        if len(I_of_q_by_frame) > 1:
            label = f"frame_{frame_number+1}"
            save_suffix = f"_frame_{frame_number+1}"
        q_min = np.min(result.mod_q)
        q_max = np.max(result.mod_q)

        if linear_binning:
            q_bins = determine_1d_linear_bins(q_min, q_max, numQBins1D)
        else:
            if (flag_logqbinsperdecade == '') or (flag_logqbinsdecadecenter is None):
                q_bins = determine_1d_log_bins(q_min, q_max, n_bins=numQBins1D,
                                               n_bins_per_decade=None,
                                               decade_on_center=flag_logqbinsdecadecenter,
                                               even_decade=flag_logqbinsevendecade)
            else:
                q_bins = determine_1d_log_bins(q_min, q_max, n_bins=None,
                                               n_bins_per_decade=flag_logqbinsperdecade,
                                               decade_on_center=flag_logqbinsdecadecenter,
                                               even_decade=flag_logqbinsevendecade)

        # [CD, 2/10/2020] added weighting option. default should be NOWEIGHT
        if flag_weighted:
            binned_i_of_q = bin_intensity_into_q1d(result, q_bins, BinningMethod.WEIGHTED)
            msapi.logger.notice('...BinningMethod = WEIGHTED.')
        else:
            binned_i_of_q = bin_intensity_into_q1d(result, q_bins, BinningMethod.NOWEIGHT)
            msapi.logger.notice('...BinningMethod = NOWEIGHT.')
        log_binned_i_of_q[frame_number] = copy.deepcopy(binned_i_of_q)

        # [CD, 1/30/2020] do we want to have an option to change btn weighted and noweighted
        # issue 322
        ax.errorbar(binned_i_of_q.mod_q, binned_i_of_q.intensity, yerr=binned_i_of_q.error, label=label)
        filename = os.path.join(json_conf["outputDir"],
                                json_params["outputFilename"] + save_suffix + '_Iq.txt')
        # save_ascii_binned_1D(filename, title, binned_i_of_q)
        # [CD, 2/10/2020] save_iqmod is preferred
        save_iqmod(binned_i_of_q, filename, float_format='%.6E')  # sep = ' ', float_format='%.6f'
    if label:
        ax.legend()
    ax.set_ylabel("Intensity")
    ax.set_xlabel("$Q (\AA^{-1})$")  # noqa W605
    ax.set_xscale('log')
    ax.set_yscale('log')

    suffix = ".png"
    picture_file = os.path.join(json_conf["outputDir"],
                                json_params["outputFilename"] + '_Iq' + suffix)
    # [CD 2/10/2020] added _Iq to the png filename
    fig.savefig(picture_file)

    # 2D
    frame_label = ''
    log_iqxqy = {}
    for frame_number, result in enumerate(I_of_qxqy_by_frame):
        if len(I_of_q_by_frame) > 1:
            frame_label = f"_frame_{frame_number+1}"
        iqxqy = get_Iqxqy(result, json_params["configuration"]["outputDir"],
                          json_params["outputFilename"], label=frame_label,
                          nbins=numQBins2D,
                          weighting=flag_weighted)
        log_iqxqy[frame_number] = iqxqy
        # [CD, 1/30/2020] option btn weighted and noweighted ?
        # [CD, 2/10/2020] option for weighting has been added

    # list of arguments for log file =======================================================
    filename = os.path.join(json_params["configuration"]["outputDir"], output_file + '_reduction_log.hdf')
    starttime = datetime.now().isoformat()
    # username = 'Neymar'
    pythonfile = __file__
    reductionparams = log_json_params
    specialparameters = {'beam_center': {'x': config['center_x'],
                                         'y': config['center_y'],
                                         },
                         'sample_transmission': sample_transmission_dict,
                         'background_transmission': background_transmission_dict,
                         }
    detectordata = {}
    for _key in log_binned_i_of_q.keys():
        name = "frame_{}".format(_key+1)
        detectordata[name] = {'iq': log_binned_i_of_q[_key],
                              'iqxqy': log_iqxqy[_key]}
    samplelogs = {'main': SampleLogs(ws)}
    drtsans.savereductionlog(filename=filename,
                             detectordata=detectordata,
                             reductionparams=reductionparams,
                             pythonfile=pythonfile,
                             starttime=starttime,
                             specialparameters=specialparameters,
                             samplelogs=samplelogs,
                             )

    # [CD 2/7/2020] log 'finish'
    msapi.logger.notice('...Reduction finished.')
