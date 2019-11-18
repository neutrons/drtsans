""" Common utility functions for all SANS """
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import drtsans  # noqa E402
import mantid.simpleapi as msapi  # noqa E402
from drtsans.iq import BinningMethod, BinningParams  # noqa E402
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402


def setup_configuration(json_params, instrument):
    """
        Extract configuration from json file passed through Shaman
    """
    # Currently unused
    # json_params['configuration']['sampleApertureSize']
    # json_params['configuration']['useDetectorTubeType']
    # json_params['configuration']['useThetaDepTransCorrection']
    # json_params['configuration']['nPixelsRadiusForTransmission']

    config = dict(is_wing=False,
                  dark_current=None,
                  sensitivity_file_path=None,
                  center_x=None,
                  center_y=None,
                  center_y_wing=None,
                  detector_offset=0,
                  sample_offset=0,
                  mask_detector=None,
                  flux_method='time',
                  solid_angle=True,
                  mask=None,
                  mask_panel=None,
                  )

    # Dark current
    if json_params['configuration']['useDarkFileName']:
        config['dark_current'] = json_params['configuration']['darkFileName']

    # Sensitivity
    if json_params['configuration']['useSensitivityFileName']:
        config['sensitivity_file_path'] = json_params['configuration']['sensitivityFileName']

    # Normalization
    config['flux_method'] = json_params['configuration']['normalization'].lower()

    # Offsets
    if json_params['configuration']['useSampleOffset']:
        config['sample_offset'] = json_params['configuration']['sampleOffset']

    if json_params['configuration']['useDetectorOffset']:
        config['detector_offset'] = json_params['configuration']['detectorOffset']

    # Solid angle
    config['solid_angle'] = json_params['configuration']['useSolidAngleCorrection']

    # Masking
    # TODO: get default mask from configuration
    if json_params['configuration']["useDefaultMask"]:
        default_mask = []
        w = msapi.LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=uwd())
        for d in default_mask:
            msapi.MaskBTP(Workspace=w, **d)
        config["mask"] = msapi.ExtractMask(w, OutputWorkspace=uwd()).OutputWorkspace
    elif json_params['configuration']['useMaskFileName']:
        config["mask"] = json_params['configuration']["maskFileName"]

    if json_params['configuration']["useMaskBackTubes"]:
        config["mask_panel"] = "back"

    return config


def get_Iq(q_data, output_dir, output_file, label='', linear_binning=True, nbins=100):
    """
        Compute I(q) from corrected workspace
    """
    q_min = np.min(q_data.mod_q)
    q_max = np.max(q_data.mod_q)
    bin_params = BinningParams(min=q_min, max=q_max, bins=nbins)
    iq_output = drtsans.iq.bin_intensity_into_q1d(q_data,
                                                  bin_params=bin_params,
                                                  linear_binning=linear_binning,
                                                  bin_method=BinningMethod.NOWEIGHT)

    filename = os.path.join(output_dir, output_file + label + '_Iq.txt')
    save_ascii_binned_1D(filename, "I(Q)", iq_output)

    fig, ax = plt.subplots()
    ax.errorbar(iq_output.mod_q, iq_output.intensity, yerr=iq_output.error, label="I(Q)")
    ax.set_xlabel("$Q (\AA^{-1})$")  # noqa W605
    plt.ylabel('Intensity')
    filename = os.path.join(output_dir, output_file + label + '_Iq.png')
    fig.savefig(filename)


def get_Iqxqy(q_data, output_dir, output_file, label='', nbins=100):
    """
        Compute I(qx,qy) from corrected workspace
    """
    qx_min = np.min(q_data.qx)
    qx_max = np.max(q_data.qx)
    binning_x = BinningParams(qx_min, qx_max, nbins)
    qy_min = np.min(q_data.qy)
    qy_max = np.max(q_data.qy)
    binning_y = BinningParams(qy_min, qy_max, nbins)

    iq_output = drtsans.iq.bin_iq_into_linear_q2d(q_data,
                                                  qx_bin_params=binning_x,
                                                  qy_bin_params=binning_y,
                                                  method=BinningMethod.NOWEIGHT)

    filename = os.path.join(output_dir, output_file + label + '_Iqxqy.txt')
    save_ascii_binned_2D(filename, "I(Qx,Qy)", iq_output)

    fig, ax = plt.subplots()
    pcm = ax.pcolormesh(iq_output.qx, iq_output.qy, iq_output.intensity.T,
                        norm=colors.LogNorm())
    fig.colorbar(pcm, ax=ax)
    picture_file = os.path.join(output_dir, output_file + label + '_Iqxqy.png')
    fig.savefig(picture_file)
