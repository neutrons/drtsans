"""
    GPSANS: preparing sensitivities script
    This implements https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/205
"""
import json
import os
import sys
from drtsans.mono.normalization import normalize_by_monitor
from drtsans.mono.load import load_events
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
import drtsans.mono.gpsans as gp
from drtsans.mono.gpsans.prepare_sensitivity import prepare_sensitivity
from drtsans.process_uncertainties import set_init_uncertainties
import numpy as np
import matplotlib.pyplot as plt
from mantid.api import AnalysisDataService as mtd
import mantid.simpleapi as msapi  # noqa E402
from mantid.simpleapi import SaveNexusProcessed
from mantid.simpleapi import CreateWorkspace, MaskDetectors
import h5py
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


DEBUGMODE = True


def export_detector_view(ws, png_name):
    """Export detector view to a PNG file

    This method is for debugging purpose

    Parameters
    ----------
    ws : ~mantid.api.MatrixWorkspace
        Workspace to plot in detector view
    png_name : str
        Path of the output PNG file

    Returns
    -------
    None

    """
    if isinstance(ws, np.ndarray):
        vec_y = ws
    else:
        vec_y = ws.extractY().transpose()
    vec_y = vec_y.reshape((192, 256)).transpose()

    plt.imshow(vec_y)

    plt.savefig(png_name)

    return


def load_data(nexus_file_name):
    """Load data

    Parameters
    ----------
    nexus_file_name : int, str
        Examples: ``55555`` or ``CG3_55555`` or file path.

    Returns
    -------
    ~mantid.dataobjects.EventWorkspace
        Normalized event workspace where data file is loaded to

    """
    ws_name = os.path.basename(nexus_file_name).split('.')[0]
    ws = load_events(nexus_file_name, output_workspace=ws_name, data_dir=None, overwrite_instrument=False)

    # Normalize by monitor
    ws = normalize_by_monitor(ws)

    return ws


def mask_data(data_ws, beam_center_ws):
    """Mask detectors in a workspace

    Mask (1) beam center, (2) top and (3) bottom

    Parameters
    ----------
    data_ws : ~mantid.api.MatrixWorkspace
        Flood data workspace
    beam_center_ws : ~mantid.api.MatrixWorkspace
        Beam center workspace used to generate beam center mask

    Returns
    -------
    ~mantid.api.MatrixWorkspace

    """
    # Use beam center ws to find beam center
    xc, yc = gp.find_beam_center(beam_center_ws)

    # Center detector to the data workspace (change in geometry)
    gp.center_detector(data_ws, xc, yc)

    # Mask the new beam center by 65 mm (Lisa's magic number)
    det = list(circular_mask_from_beam_center(data_ws, 65))
    masks = apply_mask(data_ws, mask=det)
    print('DEBUG masks: {}'.format(masks))

    # Mask top and bottom: both 8 rows
    data_ws_name = data_ws.name()
    apply_mask(data_ws, Pixel='1-8,249-256')
    data_ws = mtd[data_ws_name]

    if DEBUGMODE:
        # Optionally as a debugging step, save the workspace with masked applied
        SaveNexusProcessed(InputWorkspace=data_ws, Filename='{}_masked.nxs'.format(data_ws))
        export_detector_view(data_ws, '{}_masked.png'.format(data_ws))

    return data_ws


def setup_configuration(json_params):
    """Extract configuration

    Parameters
    ----------
    json_params : str
        Json parameter dictionary

    Returns
    -------

    """
    # Parse ITPS and runs
    ipts_number = json_params['IPTS-Number']
    flood_runs = json_params['Flood Runs']
    direct_beam_runs = json_params['Direct Beam Runs']

    # Check
    if len(flood_runs) != len(direct_beam_runs):
        raise RuntimeError('Number of flood files ({}) shall be equal to number of direct beam files ({})'
                           ''.format(len(flood_runs), len(direct_beam_runs)))

    # Form NeXus file names from IPTS and run
    flood_files = list()
    direct_beam_files = list()

    for f_index in range(len(flood_runs)):
        flood_files.append('/HFIR/CG2/IPTS-{}/nexus/CG2_{}.nxs.h5'.format(ipts_number, flood_runs[f_index]))
        direct_beam_files.append('/HFIR/CG2/IPTS-{}/nexus/CG2_{}.nxs.h5'.format(ipts_number,
                                                                                direct_beam_runs[f_index]))
    # END-FOR

    # Monitor
    monitor_array = np.array(json_params['Monitor'])

    return flood_files, direct_beam_files, monitor_array


def prepare_data(flood_run_ws_list, beam_center_run_ws_list):
    """Prepare data

    Find beam center for each flood file
    Mask top and bottom detectors
    Set uncertainties to each file
    Remove direct beam from flood

    :param flood_run_ws_list:
    :param beam_center_run_ws_list:
    :return:
    """
    num_ws_pairs = len(flood_run_ws_list)
    num_spec = flood_run_ws_list[0].getNumberHistograms()

    masked_flood_list = [None] * num_ws_pairs
    for i_pair in range(num_ws_pairs):
        # use beam center run to locate the beam center and then do mask to data workspace
        bc_ws = beam_center_run_ws_list[i_pair]
        flood_ws = flood_run_ws_list[i_pair]
        # mask
        masked_flood_ws = mask_data(flood_ws, bc_ws)
        # set uncertainties:
        # output: masked are zero intensity and zero error
        masked_flood_ws = set_init_uncertainties(flood_ws, masked_flood_ws)
        # append
        masked_flood_list[i_pair] = masked_flood_ws
    # END-FOR

    # Combine to numpy arrays: N, M
    flood_array = np.ndarray(shape=(num_ws_pairs, num_spec), dtype=float)
    sigma_array = np.ndarray(shape=(num_ws_pairs, num_spec), dtype=float)
    for f_index in range(num_ws_pairs):
        flood_array[f_index][:] = masked_flood_list[f_index].extractY().transpose()[0]
        sigma_array[f_index][:] = masked_flood_list[f_index].extractE().transpose()[0]

    # Convert all to NaN
    masked_items = np.where(sigma_array < 1E-16)

    # set values
    flood_array[masked_items] = np.nan
    sigma_array[masked_items] = np.nan

    return flood_array, sigma_array


def save_to_nexus(parent_ws, sensitivities, sensitivities_error, output_file):
    """

    Parameters
    ----------
    parent_ws : ~mantid.api.MatrixWorkspace
        Parent GSPANS workspace serving as
    sensitivities
    sensitivities_error
    output_file

    Returns
    -------

    """
    # Create a workspace for sensitivities
    vec_x = parent_ws.extractX().flatten()
    num_spec = parent_ws.getNumberHistograms()
    sens_ws_name = 'sensitivities'

    nexus_ws = CreateWorkspace(DataX=vec_x, DataY=sensitivities, DataE=sensitivities_error,
                               NSpec=num_spec, ParentWorkspace=parent_ws,
                               OutputWorkspace=sens_ws_name)

    # Check masks
    num_masked_pixels = 0
    for i in range(nexus_ws.getNumberHistograms()):
        if nexus_ws.getDetector(i).isMasked():
            num_masked_pixels += 1
    print('Original masked: {}'.format(num_masked_pixels))

    # Mask
    masked_ws_list = np.where(np.isnan(sensitivities))[0]
    MaskDetectors(Workspace=sens_ws_name, WorkspaceIndexList=masked_ws_list)

    # Check masks
    num_masked_pixels = 0
    for i in range(nexus_ws.getNumberHistograms()):
        if nexus_ws.getDetector(i).isMasked():
            num_masked_pixels += 1
    print('Output masked: {}'.format(num_masked_pixels))

    # Export
    SaveNexusProcessed(InputWorkspace=sens_ws_name, Filename=output_file)

    return


def main(argv):
    """Main function

    Parameters
    ----------
    argv

    Returns
    -------

    """
    if len(argv) < 2:
        raise RuntimeError("Prepare sensitivities requires a parameter json string")

    # Parse the input JSON string for configuration
    if os.path.isfile(argv[1]):
        # case: Input is a Json file
        print(argv[1])
        with open(argv[1], 'r') as fd:
            json_params = json.load(fd)
    else:
        # case: Input is a Json string
        json_string = " ".join(argv[1:])
        json_params = json.loads(json_string)
    msapi.logger.notice(json.dumps(json_params, indent=2))

    output_file = json_params['output_file_name']

    # Set up the configuration
    flood_files, direct_beam_files, monitor_array = setup_configuration(json_params)

    # Load data/runs
    flood_ws_list = list()
    beam_center_ws_list = list()
    num_ws_pairs = len(flood_files)
    for f_index in range(num_ws_pairs):
        # load flood
        flood_ws_i = load_data(flood_files[f_index])
        flood_ws_list.append(flood_ws_i)

        # load beam center
        beam_center_ws_i = load_data(direct_beam_files[f_index])
        beam_center_ws_list.append(beam_center_ws_i)
    # END-FOR

    flood_data_array, flood_sigma_array = prepare_data(flood_ws_list, beam_center_ws_list)
    # Export flood data and sigma
    sens_h5 = h5py.File('{}.h5'.format(output_file.split('.')[0]), 'w')
    flood_group = sens_h5.create_group('Flood')
    flood_group.create_dataset('flood', data=flood_data_array)
    flood_group.create_dataset('flood error', data=flood_sigma_array)

    # Calculate sensitivities for each file
    sens_set = prepare_sensitivity(flood_data_matrix=flood_data_array, flood_sigma_matrix=flood_sigma_array,
                                   monitor_counts=monitor_array,
                                   threshold_min=0.5, threshold_max=1.5)

    # Export sensitivities calculated in file to ....
    sensitivities, sensitivities_error = sens_set

    # Export 2D view
    export_detector_view(sensitivities, 'Sensitivity.png')
    # Export to hdf5
    # sens_h5 = h5py.File('{}.h5'.format(output_file.split('.')[0]), 'w')
    sens_group = sens_h5.create_group('Sensitivities')
    sens_group.create_dataset('sensitivities', data=sensitivities)
    sens_group.create_dataset('sensitivities error', data=sensitivities_error)
    sens_h5.close()

    # Export to the NeXus file
    save_to_nexus(beam_center_ws_list[0], sensitivities, sensitivities_error, output_file)

    return


def generate_test_json():
    """Generate test case

    From the commission runs in Nov. 2019 provided by instrument scientist

    585 Flood 6m center 2019/10/31 20:04:31 EDT 2:00:00 6.08 5.36 8
    586 db for 6m center 2019/10/31 22:06:11 EDT 0:02:00 6.08 5.36 8
    587 Flood 6 m 200 off center 2019/10/31 22:12:28 EDT 2:00:00 6.08 5.36 8
    588 db 6 m 200 off center 2019/11/01 00:14:08 EDT 0:02:00 6.08 5.36 8
    589 Flood 6 m 400 off center 2019/11/01 00:20:27 EDT 2:00:09 6.08 5.36 8
    590 db 6 m 400 off center 2019/11/01 02:22:07 EDT 0:02:00 6.08 5.36 8

    Returns
    -------
    str
        A standard JSON string (dict)

    """
    # Set up the dictionary for JSON

    json_dict = {'IPTS-Number': 24664,
                 'Flood Runs': [1697, 1701, 1699],
                 'Direct Beam Runs': [1698, 1702, 1700],
                 'Monitor': [1., 1., 1],
                 'output_file_name': 'test_gp_sensitivity.nxs'}
    json_str = json.dumps(json_dict)

    return json_str


if __name__ == "__main__":
    # set up a test case if there is no JSON filled in
    if len(sys.argv) == 1:
        script_argv = ['', generate_test_json()]
    else:
        script_argv = sys.argv

    # Call main
    main(script_argv)
