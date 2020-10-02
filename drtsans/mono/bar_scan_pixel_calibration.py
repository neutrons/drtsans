from copy import deepcopy
from matplotlib import pyplot as plt
import numpy as np
import os
import time
# Mantid imports
from mantid.api import mtd
from mantid.simpleapi import (CreateWorkspace, DeleteWorkspaces, LoadEventNexus, LoadNexus,
                              HFIRSANS2Wavelength, SaveNexus)
# drtsans imports
from drtsans.mono.gpsans import (apply_calibrations, apply_mask, calculate_apparent_tube_width,
                                 calculate_barscan_calibration)
from drtsans.pixel_calibration import Table
from drtsans.tubecollection import TubeCollection
import sys


def locate_bar_scan_files(ipts_number, exp_number, scan_number, pt_numbers, root_dir):
    # Convert from SPICE xml to event Nexus
    ipts_directory = os.path.join(root_dir, f'IPTS-{ipts_number}/shared/Exp{exp_number}/')
    if os.path.exists(ipts_directory) is False:
        os.mkdir(ipts_directory)

    data_files = list()
    err_msg = ''
    for pt_number in pt_numbers:  # bar_scan_files.items():
        event_nexus_name = 'CG2_{:04}{:04}{:04}.nxs.h5'.format(exp_number, scan_number, pt_number)
        event_nexus_name = os.path.join(ipts_directory, event_nexus_name)
        if not os.path.exists(event_nexus_name):
            err_msg += f'Pt {pt_number} has been not been converted to {event_nexus_name} yet\n'
        else:
            data_files.append(event_nexus_name)
    # END-FOR
    if len(err_msg) > 0:
        print(f'[ERROR] {err_msg}')
        sys.exit(-1)

    return data_files


def locate_flood_file(flood_ipts_number, flood_exp_number, flood_scan_number, flood_pt_number, root_dir):
    # convert flood file
    # TODO: refactor convert XML to event Nexus for all files
    flood_ipts_directory = f'/HFIR/CG2/IPTS-{flood_ipts_number}/shared/Exp{flood_exp_number}/'
    if os.path.exists(flood_ipts_directory) is False:
        os.mkdir(flood_ipts_directory)
    flood_nexus_name = 'CG2_{:04}{:04}{:04}.nxs.h5'.format(flood_exp_number, flood_scan_number, flood_pt_number)
    flood_nexus_name = os.path.join(flood_ipts_directory, flood_nexus_name)
    if not os.path.exists(flood_nexus_name):
        sys.exit(-1)

    return flood_nexus_name


def generate_intermediate_files(bar_scan_files, save_dir):
    print('####\n\nCreating intermediate files, one for each barscan run. This can take up to one hour')

    barscan_dataset = list()  # list of histogram files
    # Convert to histogram for future use with efficiency
    for bar_scan_run in bar_scan_files:
        base_name = os.path.basename(bar_scan_run).split('.')[0]
        file_histogram = os.path.join(save_dir, f'{base_name}.nxs')
        if os.path.exists(file_histogram):
            # exist
            print('File {} already exists'.format(file_histogram))
        else:
            # load and recreate
            workspace_events = f'{base_name}_events'
            LoadEventNexus(Filename=bar_scan_run, LoadMonitors=False, OutputWorkspace=workspace_events)

            workspace_counts = f'{base_name}_counts'
            HFIRSANS2Wavelength(InputWorkspace=workspace_events, OutputWorkspace=workspace_counts)

            SaveNexus(InputWorkspace=workspace_counts, Filename=file_histogram)
            # Clean workspace
            DeleteWorkspaces([workspace_events, workspace_counts])

        # append file name for next step
        barscan_dataset.append(file_histogram)

    return barscan_dataset


def generate_pixel_map(bar_scan_files, flood_file, save_dir_root, database_file_base,
                       mask_file='testdata/mask_pixel_map.nxs'):
    """Generate pixel calibration map from bar scan

    Returns
    -------
    generator
        bar scan data set, calibration (step 1),

    """
    # Check output directory
    if not os.path.exists(save_dir_root):
        os.mkdir(save_dir_root)
    save_dir = os.path.join(save_dir_root, f'{len(bar_scan_files)}_runs')
    print(f'save to {save_dir}')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir, exist_ok=True)

    # Create intermediate files
    barscan_dataset = generate_intermediate_files(bar_scan_files, save_dir)
    # return bar scan data set
    yield barscan_dataset

    print('#####\n\nCalculating the barscan calibration with the default formula. This takes ~10 minutes')
    start_time = time.time()
    # TODO - THIS IS SO MAGIC!
    formula = '565 - {y} - 0.0914267 * (191 - {tube})'
    calibration = calculate_barscan_calibration(barscan_dataset, formula=formula)
    print('Calibration took ', int((time.time() - start_time) / 60), 'minutes')

    print('####\n\nRemoving the Bar Tilt and Centering the Detector')
    calibration = untilt_and_center(calibration)
    # return first calibration
    calibration.state_flag = 1
    yield calibration  # calibration step 1

    print('#####\n\nSaving the calibration .. May overwrite already saved calibration')
    # Notice we overwrite the already saved calibration, which will happen if we run this notebook more than once.

    # table file name (find out)
    database_file = os.path.join(save_dir, database_file_base)
    print(f'Database file {database_file}')
    cal_dir = os.path.join(os.path.dirname(database_file), 'tables')  # directory where to save the table file
    os.makedirs(cal_dir, exist_ok=True)  # Create directory, and don't complain if already exists
    test_table_file = os.path.join(cal_dir, 'Test_CG2_Pixel_Calibration' + '.nxs')
    # save
    calibration.save(overwrite=True, database=database_file, tablefile=test_table_file)

    print('#####\n\napply the calibration to the flood run as a test')
    flood_ws_name = 'flood_run'
    LoadEventNexus(Filename=flood_file, OutputWorkspace=flood_ws_name)
    HFIRSANS2Wavelength(InputWorkspace=flood_ws_name, OutputWorkspace=flood_ws_name)
    if False:
        calibrated_flood_ws_name = f'{flood_ws_name}_calibrated'
        apply_calibrations(flood_ws_name, calibrations='BARSCAN', output_workspace=calibrated_flood_ws_name,
                           database=database_file)

    print('#####\n\nCalculating the Tube Width Calibration')

    # Mask file containing the detector ID's comprising the beam center.
    # mask_file = f'/HFIR/CG2/IPTS-{ipts}/shared/pixel_flood_mask.nxs'

    apply_mask(flood_ws_name, mask=mask_file)
    start_time = time.time()
    calibration = calculate_apparent_tube_width(flood_ws_name, load_barscan_calibration=True,
                                                db_file=database_file)
    print('Calibration took ', int(time.time() - start_time), 'seconds')

    print('#####\n\nSaving the Tube Width calibration')
    # Notice we overwrite the already saved calibration, which will happen if we run this notebook more than once.
    # calibration.save(overwrite=True)
    calibration.save(overwrite=True, database=database_file, tablefile=test_table_file)
    calibration.state_flag = 2
    yield calibration, flood_ws_name  # calibration 2

    yield test_table_file


def generate_spice_pixel_map(ipts_number, exp_number, scan_number, pt_numbers,
                             flood_ipts_number, flood_exp_number, flood_scan_number, flood_pt_number,
                             root_dir, save_dir_root):
    bar_scan_files = locate_bar_scan_files(ipts_number, exp_number, scan_number, pt_numbers, root_dir)
    flood_file = locate_flood_file(flood_ipts_number, flood_exp_number, flood_scan_number, flood_pt_number, root_dir)

    # data base file name
    database_file_dict = {'CG2': '/HFIR/CG2/shared/calibration/pixel_calibration.json',
                          'CG3': '/HFIR/CG3/shared/calibration/pixel_calibration.json'}

    database_file_base = os.path.basename(database_file_dict['CG2'])

    return generate_pixel_map(bar_scan_files, flood_file, save_dir_root, database_file_base)


def plot_histograms(input_workspace, legend=[], xlabel='X-axis', ylabel='Y-axis', title='', linewidths=[]):
    r"""Line plot for the histograms of a workspace"""
    workspace = mtd[str(input_workspace)]
    number_histograms = workspace.getNumberHistograms()
    if len(legend) != number_histograms:
        legend = [str(i) for i in range(number_histograms)]
    if len(linewidths) != number_histograms:
        linewidths = [1] * number_histograms
    fig, ax = plt.subplots(subplot_kw={'projection': 'mantid'})
    for workspace_index in range(number_histograms):
        ax.plot(workspace, wkspIndex=workspace_index, label=legend[workspace_index],
                linewidth=linewidths[workspace_index])
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='out')
    ax.grid(True)
    fig.show()


def linear_density(workspace):
    r"""
    Tube total intensity per non-masked pixel and per unit length of tube width

    Integrate the total intensity per tube and divide by the number of non-masked pixels in the tube,
    and by the tube width. Front end tubes collect more intentity than the back tubes.
    Similarly, front end tubes have a larger apparent tube width than back tubes.
    The ratio of total intensity to width should be similar for front and end tubes after the calibration.
    """
    collection = TubeCollection(workspace, 'detector1').sorted(view='fbfb')
    intensities = np.array([np.sum(tube.readY) for tube in collection])
    widths = np.array([tube.width for tube in collection])
    number_pixels_not_masked = np.array([np.sum(~tube.isMasked) for tube in collection])
    return list(intensities / (number_pixels_not_masked * widths))


def report_tilt(pixel_positions):
    r"""
    Variation in the position of the top and bottom pixels as a function of tube index.
    We perform a linear regression of this variation.
    """
    # Create a 2D array of pixel heights, dimensions are (number_tubes x pixels_in_tube)
    pixel_in_tube_count = 256
    tube_count = int(len(pixel_positions) / pixel_in_tube_count)
    positions = np.array(pixel_positions).reshape((tube_count, pixel_in_tube_count))

    def fit(tube_tip_positions):
        r"""This function will fit the bottom or top pixels against the tube index"""
        tube_indexes = np.arange(tube_count)  # heights as function of tube index
        coeffs = np.polyfit(tube_indexes, tube_tip_positions, 1)
        fitted = np.poly1d(coeffs)(tube_indexes)  # fitted positions of the tube tip
        return coeffs, fitted

    for location, tip_positions in (['top', positions[:, -1]], ['bottom', positions[:, 0]]):
        coeffs, fitted = fit(tip_positions)  # fit against tube index
        # Plot the raw positions and the fitted positions
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(np.arange(tube_count), tip_positions)
        ax.plot(np.arange(tube_count), fitted)
        ax.set_title(f'{location} pixels')
        # Print a few representative properties of the tilt
        print(location, ' pixels:')
        print(f'    slope = {1000 * coeffs[0]:.3f} mili-meters / tube')
        print(
            f'    position difference between last and first tube = {1000 * (fitted[-1] - fitted[0]):.3f} mili-meters')
        # save figure
        plt.savefig(f'{location}.png')


def untilt_and_center(a_calibration):
    r"""
    Removing the Bar Tilt and Centering the Detector

    Thinking of the fitted positions for the bottom and top pixels, we can think of the detector array
    as a deformed rectangle (angles between sides different than 90 degrees), which must be transformed
    into a rectangle with squared angles (angles between sides equal to 90 degrees).

    We take the tube in the middle of the main detector array as our reference. We will adjust
    every other tube so that for every tube, its top and bottom *fitted* pixel positions
    will coincide with the top and bottom *fitted* positions of the middle tube.

    Also, since top and bottom fitted positions have a different variation with tube index,
    the fitted tube lenght changes sligtly with tube index. Thus, we will rescale the fitted
    tube length to coincide with the fitted tube length of the middle tube. This amounts to
    a rescaling of pixel heights.

    Finally, after removing the tilt we displace the detector so that the center of mass lies at `Y=0`.
    """
    # Create a 2D array of pixel heights, dimensions are (number_tubes x pixels_in_tube)
    pixel_in_tube_count = 256
    tube_count = int(len(a_calibration.positions) / pixel_in_tube_count)
    positions = np.array(a_calibration.positions).reshape((tube_count, pixel_in_tube_count))
    heights = np.array(a_calibration.heights).reshape((tube_count, pixel_in_tube_count))

    def fit(tube_tip_positions):
        r"""This function will fit the bottom or top pixels against the tube index"""
        tube_indexes = np.arange(tube_count)  # heights as function of tube index
        coeffs = np.polyfit(tube_indexes, tube_tip_positions, 1)
        fitted = np.poly1d(coeffs)(tube_indexes)  # fitted positions of the tube tip
        return coeffs, fitted

    _, fitted_top = fit(positions[:, -1])  # fitted positions of the tube tops
    _, fitted_bottom = fit(positions[:, 0])  # fitted positions of the tube bottom
    # We'll adjust the positions of the tubes to comply with the middle tube
    tube_reference_index = int(tube_count / 2)  # tube in the middle of the detector
    tube_length_reference = fitted_top[tube_reference_index] - fitted_bottom[tube_reference_index]
    # shifts_top indicate the difference in fitted positions for the tube tops with respect to the fitted positions
    # for the top of the middle tube
    shifts_top = fitted_top[tube_reference_index] - fitted_top
    shifts_bottom = fitted_bottom[tube_reference_index] - fitted_bottom
    # Calculate now the shifts for every single pixel, going tube by tube
    pixel_indexes = np.arange(pixel_in_tube_count)
    shifts = list()
    scalings = list()
    for tube_index in range(tube_count):
        a, b = shifts_bottom[tube_index], shifts_top[tube_index]
        shifts_in_tube = a + (b - a) * pixel_indexes / pixel_in_tube_count
        shifts.append(shifts_in_tube)
        tube_length = fitted_top[tube_index] - fitted_bottom[tube_index]
        scalings_in_tube = [tube_length_reference / tube_length] * pixel_in_tube_count
        scalings.append(scalings_in_tube)

    positions_new = positions + np.array(shifts)
    heights_new = heights * np.array(scalings)

    # Set CM at y=0
    positions_new -= np.mean(positions.ravel())

    # retrieve components from the main calibration in order to construct a new calibration
    metadata = deepcopy(a_calibration.metadata)
    detector_ids = deepcopy(a_calibration.detector_ids)
    recalibration = Table(metadata,
                          detector_ids=detector_ids,
                          positions=positions_new.ravel(),
                          heights=heights_new.ravel())
    return recalibration

