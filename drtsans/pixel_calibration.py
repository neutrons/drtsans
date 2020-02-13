import numpy as np
import numexpr
import tinydb
import collections

r""" Hyperlinks to mantid algorithms
CloneWorkspace <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
DeleteWorkspaces <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspaces-v1.html>
Integration <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
LoadInstrument <https://docs.mantidproject.org/nightly/algorithms/LoadInstrument-v1.html>
MaskDetectors <https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html>
MaskDetectorsIf <https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html>
ReplaceSpecialValues <https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html>
"""
from mantid.simpleapi import (CloneWorkspace, CreateWorkspace, DeleteWorkspaces, FilterEvents,
                              GenerateEventsFilter, Integration, Load, LoadInstrument, MaskDetectors,
                              MaskDetectorsIf, ReplaceSpecialValues)
from mantid.api import mtd

r"""
Hyperlinks to drtsans functions
namedtuplefy, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
SampleLogs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.instruments import InstrumentEnumName, instrument_enum_name, instrument_standard_name
from drtsans.path import exists as file_exists
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.tubecollection import TubeCollection


__all__ = ['calculate_pixel_calibration', 'load_and_apply_pixel_calibration']

# flags a problem identifying the pixel corresponding to the bottom of the shadow cast by the bar
INCORRECT_PIXEL_ASSIGNMENT = -1


def _consecutive_true_values(values, how_many, reverse=False, raise_message=None):
    r"""
    Find first array index of consecutive `how_many` True values.

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    values: list
        list of `True` and `False` items
    how_many: int
        Number of desired consecutive `True` values
    raise_message: str
        Exception message. No exception if :py:obj:`None`, but INCORRECT_PIXEL_ASSIGNMENT is returned

    Returns
    -------
    int

    Raises
    ------
    IndexError
        If no index is found for any of the edges
    RuntimeError
        If a faulty tube is found
    """
    # use the array or the reverse one
    truth_array = values[::-1] if reverse else values
    # create a sub-array of length how_many of True values that we want to find
    pattern = [True]*how_many
    # loop over the input data and return the first index where the next
    # how_many elements match the pattern
    for i in range(len(truth_array) - how_many):
        if truth_array[i: i + how_many] == pattern:
            return len(values) - i - 1 if reverse else i
    # raise an error if the pattern is not found
    else:
        if raise_message is not None:
            raise IndexError(raise_message)
        return INCORRECT_PIXEL_ASSIGNMENT  # signal for non-identified value


@namedtuplefy
def find_edges(intensities, tube_threshold=0.2, shadow_threshold=0.3,
               tube_edge_min_width=3, shadow_edge_min_width=4,
               min_illuminated_length=7):
    r"""
    Find the active length of the tube and the shadow region

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    intensities: list
        pixel pixel_intensities along the tube.
    tube_threshold: float
        fraction of the average intensity to determine the tube edges.
    shadow_threshold: float
        fraction of the average intensity to determine the shadow region.
    tube_edge_min_width: int
        required minimum number of consecutive good pixels above
        the tube threshold
    shadow_edge_min_width: int
        required minimum number of consecutive shadowed pixels
    min_illuminated_length: int
        minimum number of illuminated pixels on the active length

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - bottom_pixel: first illuminated pixel
        - top_pixel: last illuminated pixel
        - bottom_shadow_pixel: first shadowed pixel
        - above_shadow_pixel= first illuminated pixel above the shadow region
    """
    # calculate minimum intensity thresholds for tube ends and shadows
    average_intensity = np.average(intensities)
    end_threshold = tube_threshold * average_intensity
    shadow_threshold = shadow_threshold * average_intensity

    # Find edges of the tube: want at least tube_edge_min_width pixels
    # (starting from the top or bottom of a tube) that have pixel_intensities greater
    # than the threshold

    illuminated = [bool(i > end_threshold) for i in intensities]
    # The bottom pixel is the first illuminated pixel. It is required that the next tube_edge_min_width pixels are
    # also illuminated
    bottom_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                            raise_message='Could not find bottom tube edge')
    # The top pixel is the last illuminated pixel. It is required that the previous tube_edge_min_width pixels are
    # also illuminated.
    top_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                         raise_message='Could not find top tube edge', reverse=True)

    # Find the shadow region: similar to tube edges, but in this case
    # we want shadow_edge_min_width pixel_intensities less than the shadow threshold,
    # followed by at least one intensity greater than the threshold

    # The bottom pixel shadowed by the bar is the first pixel below the intensity threshold. We require that the
    # next shadow_edge_min_width are also shadowed.
    shadowed = [bool(i < shadow_threshold) for i in intensities[bottom_pixel:]]
    bottom_shadow_pixel = bottom_pixel +\
        _consecutive_true_values(shadowed, shadow_edge_min_width, raise_message='Could not find bottom shadow edge')

    # Find the first illuminated pixel above the bar.
    illuminated = [bool(i > shadow_threshold) for i in
                   intensities[bottom_shadow_pixel + shadow_edge_min_width: top_pixel + 1]]
    # Don't raise if the pixel is not found
    above_shadow_pixel = bottom_shadow_pixel + shadow_edge_min_width +\
        _consecutive_true_values(illuminated, 1, raise_message=None)

    # Check for a faulty tube: we want a certain number of pixels not in the bar shaddow
    active_tube_length = top_pixel - bottom_pixel + 1
    shadow_length = above_shadow_pixel - bottom_shadow_pixel
    if active_tube_length < min_illuminated_length + shadow_length:
        raise RuntimeError('Faulty tube found')

    return dict(bottom_pixel=bottom_pixel, top_pixel=top_pixel,
                bottom_shadow_pixel=bottom_shadow_pixel,
                above_shadow_pixel=above_shadow_pixel)


@namedtuplefy
def fit_positions(edge_pixels, bar_positions, tube_pixels=256, order=5, ignore_value=INCORRECT_PIXEL_ASSIGNMENT):
    r"""
    Fit the position and heights of the pixels in a tube. The bar_positions as a function of
    edge pixels are fitted to a nth order polynomial (by default n=5). The positions of the pixels along the
    tube are the values of the polynomial at integer points, while the heights are the derivatives.

    Description from the master requirements document, section A2.1

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    Uses :ref:`~numpy.polynomial.polynomial.polyfit`.

    devs - Andrei Savici <saviciat@ornl.gov>,

    Parameters
    ----------
    edge_pixels: list (or numpy array)
        the bottom pixel for each bar position, as found in `find_edges` function
    bar_positions: list (or numpy array)
        the bar position from the logs for each file in the bar scan
    tube_pixels: integer
        number of pixels for which to calculate positions and heights
    order: integer
        the order of polynomial to be used in the fit (default 5)
    ignore_value: int
        certain positions of the bar (close to the top and bottom of the tube) results in incorrect assignment of the
        edge pixel. In those cases it is expected that the edge pixel has a particular value that flags incorrect
        assignment. The default value is INCORRECT_PIXEL_ASSIGNMENT. These edge pixels will be
        ignored when carrying out the fit.

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - calculated_positions: calculated positions of the pixels
        - calculated_heights: calculated pixel heights
    """
    message_len = 'The positions of the bar and edge pixels have to be the same length'
    assert len(edge_pixels) == len(bar_positions), message_len

    # Ignore the incorrectly assigned edge pixels
    edge_pixels = np.array(edge_pixels)
    bar_positions = np.array(bar_positions)
    valid_edge_pixels = edge_pixels[np.where(edge_pixels != ignore_value)]
    valid_bar_positions = bar_positions[np.where(edge_pixels != ignore_value)]

    try:
        # fit the bar positions to a 5th degree polynomial in edge_pixels
        coefficients = np.polynomial.polynomial.polyfit(valid_edge_pixels, valid_bar_positions, int(order))
        # calculate the coefficients of the derivative
        deriv_coefficients = np.polynomial.polynomial.polyder(coefficients)
        # evalutae the positions
        calculated_positions = np.polynomial.polynomial.polyval(np.arange(tube_pixels), coefficients)
        # evaluate the heights
        calculated_heights = np.polynomial.polynomial.polyval(np.arange(tube_pixels), deriv_coefficients)
    except Exception:
        calculated_positions = np.ones(tube_pixels) * np.nan
        calculated_heights = np.ones(tube_pixels) * np.nan

    return dict(calculated_positions=calculated_positions, calculated_heights=calculated_heights,
                coefficients=coefficients)


# Files storing the all pixel calibrations for each instrument. Used in `load_calibration` and `save_calibration`
database_file = {InstrumentEnumName.BIOSANS: '/HFIR/CG3/shared/pixel_calibration.json',
                 InstrumentEnumName.EQSANS: '/SNS/EQSANS/shared/pixel_calibration.json',
                 InstrumentEnumName.GPSANS: '/HFIR/CG2/shared/pixel_calibration.json'}


def load_calibration(instrument, run=None, component='detector1', database=None):
    r"""
    Load pixel calibration from the database.

    Parameters
    ----------
    instrument:str
        Name of the instrument
    run: int
        Run number to resolve which calibration to use. If :py:obj:`None`, then the lates calibration will be used.
    component: str
        Name of the double panel detector array for which the calibration was performed
    database: str
        Path to database. If :py:obj:`None`, the default database is used.

    devs - Jose Borreguero <borreguerojm@ornl.gov>,

    Returns
    -------
    dict
        Dictionary with the following entries:
        - instrument, str, Name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number associated to the calibration.
        - unit: str, the units for the positions, heights, and widths. Usually to 'mm' for mili-meters.
        - positions, list, List of Y-coordinate for each pixel.
        - heights, list, List of pixel heights.
        - widths, list, A two-item list containing the apparent widths for the front and back tubes.
    """
    # Find entries with given instrument and component.
    enum_instrument = instrument_enum_name(instrument)
    if database is None:
        database = database_file[enum_instrument]
    database_server = tinydb.TinyDB(database)
    calibrations = tinydb.Query()

    def find_matching_calibration(calibration_type_query):
        matches = database_server.search(calibration_type_query &
                                         (calibrations.instrument == enum_instrument.name) &
                                         (calibrations.component == component))
        if len(matches) == 0:
            return {}  # no calibration is found
        # Sort by decreasing run number and find appropriate calibration
        matches = sorted(matches, key=lambda d: d['run'], reverse=True)
        if run is None:
            return matches[0]  # return latest calibration
        else:
            for i, match in enumerate(matches):
                if run > match['run']:
                    return matches[i-1]  # return first calibration with a smaller run number than the query run number
                elif run == match['run']:
                    return matches[i]  # strange corner case

    # Find matching barscan calibration (pixel positions and heights)
    calibration = find_matching_calibration(calibrations.heights)
    # Find matching flat-field calibration (pixel widths) and update the calibration dictionary
    calibration['widths'] = find_matching_calibration(calibrations.widths).get('widths', None)
    database_server.close()

    return calibration


def save_calibration(calibration, database=None):
    r"""
    Save a calibration to the pixel-calibrations database.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    calibration: dict
        Dictionary containing the following required keys:
        - instrument, str, Name of the instrument
        - component, str, name of the double detector array, usually "detector1"
        - run, int, run number associated to this calibration.
        - unit: str, the units for the updated pixel dimensions. Usually 'mm' for mili-meters.
        The dictionary will contain also entries for pixel positions and heights, as well as front and back
        tube widths. One can pass a dictionary containing all three entries, or a dictionary containing pixel
        positions and heights (the results of a barscan calibration) or a dictionary containing tube widths (the
        result of a flat-field calibration)
    database: str
        Path to database file. If :py:obj:`None`, the default database is used.
    """
    # Validate calibration. Check is a mapping and check for required keys
    if isinstance(calibration, collections.Mapping) is False:
        raise ValueError('The input is not a mapping object, such as a dictionary')
    if {'instrument', 'component', 'run'} in set(calibration.keys()) is False:
        raise KeyError('One or more mandatory keys missing ("instrument", "component", "run"')

    # Make sure the instrument name is the standard instrument name. For instance, "GPSANS" instead of "CG2"
    enum_instrument = instrument_enum_name(calibration['instrument'])
    calibration['instrument'] = str(enum_instrument)  # overwrite with the standard instrument name

    # Save entry in the database
    if database is None:
        database = database_file[enum_instrument]
    database_server = tinydb.TinyDB(database)
    database_server.insert(dict(calibration))
    database_server.close()


def event_splitter(barscan_file, split_workspace=None, info_workspace=None, bar_position_log='dcal_Readback'):
    r"""


    Parameters
    ----------
    barscan_file: str
        Path to barscan run file containing multiple positions of the bar.
    split_workspace: str
        Name of the table workspace to be used as event splitter workpsce in algorithm FilterEvents. If
        :py:obj:`None`, a random name will be provided.
    info_workspace: str
        Name of the table workspace to be used along ``split_workspace`` in algorithm FilterEvents. If
        :py:obj:`None`, a random name will be provided.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.
    Returns
    -------
    list
        List of bar positions
    """
    if split_workspace is None:
        split_workspace = unique_workspace_dundername()
    if info_workspace is None:
        info_workspace = unique_workspace_dundername()

    barscans_workspace = unique_workspace_dundername()
    Load(barscan_file, OutputWorkspace=barscans_workspace)

    # Find the amount by which the position of bar is shifted
    bar_positions = SampleLogs(barscans_workspace)[bar_position_log].value
    bar_delta_positions = bar_positions[1:] - bar_positions[:-1]  # list of shifts in the position of the bar
    bar_delta_positions = bar_delta_positions[bar_delta_positions > 0]  # only shifts where the bar position increases
    # the most likely shift of the bar position, within two significant figures.
    bar_step = float(np.bincount(np.round(100 * bar_delta_positions).astype('int')).argmax()) / 100.

    GenerateEventsFilter(InputWorkspace=barscans_workspace, OutputWorkspace=split_workspace,
                         InformationWorkspace=info_workspace, UnitOfTime='Nanoseconds',
                         LogName=bar_position_log, LogValueInterval=bar_step, LogValueTolerance=bar_step / 2,
                         MinimumLogValue=min(bar_positions), MaximumLogValue=max(bar_positions))

    # Read bar positions from info_workspace using the second column. That column has entries as strings of the form:
    # "Log.dcal_Readback.From.{min}.To.{max}.Value-change-direction:both"
    bar_positions = list()
    for min_max_bar_position in mtd[info_workspace].column(1):
        min_bar_position = float(min_max_bar_position.split('.From.')[-1].split('.To.')[0])
        max_bar_position = float(min_max_bar_position.split('.To.')[-1].split('.Value')[0])
        bar_positions.append((min_bar_position + max_bar_position) / 2)

    return bar_positions


def barscan_workspace_generator(barscan_files, bar_position_log='dcal_Readback'):
    r"""
    A generator to be used in an iteration over the runs that held the bar at a fixed position, returning at each
    iteration the name of the workspace containing a single run (the set of runs make up the bar scan).

    Parameters
    ----------
     barscan_files: str, list
        Path(s) to barscan run file(s). If only one file, it should contain multiple positions of the bar.
        If a list of files, then each file contains the pixel_intensities recorded with a constant position for the
        bar.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.

    Returns
    -------
    tuple
        A two-item tuple containing, in this order:
        - the position of the bar as stated in the logs the run currently being returned.
        - the name of the workspace containing the run currently being returned.
    """
    temporary_workspaces = list()  # store the names of the workspaces to be removed

    def temporary_workspace():
        r"""Random workspace name, and flag it for removal"""
        name = unique_workspace_dundername()
        temporary_workspaces.append(name)
        return name

    if isinstance(barscan_files, str):
        spliter_workspace = temporary_workspace()
        info_workspace = temporary_workspace()
        # Table to split the barscan run into subruns, each with a fixed position of the bar
        bar_positions = event_splitter(barscan_files, split_workspace=spliter_workspace,
                                       info_workspace=info_workspace, bar_position_log=bar_position_log)
        barscans_workspace = temporary_workspace()
        Load(barscan_files, OutputWorkspace=barscans_workspace)
        splitted_workspace_group = unique_workspace_dundername()
        FilterEvents(InputWorkspace=barscans_workspace,
                     SplitterWorkspace=spliter_workspace, InformationWorkspace=info_workspace,
                     OutputWorkspaceBaseName=splitted_workspace_group,  GroupWorkspaces=True,
                     TimeSeriesPropertyLogs=[bar_position_log], ExcludeSpecifiedLogs=False)
        temporary_workspaces.append(splitted_workspace_group)
        temporary_workspaces.append('TOFCorrectWS')  # spurious workspace spawned by FilterEvents
        for i, bar_position in enumerate(bar_positions):
            yield bar_position, splitted_workspace_group + '_' + str(i)
    else:
        barscan_workspace = temporary_workspace()
        for file_name in barscan_files:
            Load(file_name, OutputWorkspace=barscan_workspace)
            bar_position = SampleLogs(barscan_workspace).find_log_with_units(bar_position_log, 'mm')
            yield bar_position, barscan_workspace
    DeleteWorkspaces(temporary_workspaces)  # clean up the now useless workspaces


def calculate_barscan_calibration(barscan_files, component='detector1', bar_position_log='dcal_Readback',
                                  formula='565 - {y}', order=5):
    r"""
    Calculate pixel positions (only Y-coordinae) as well as pixel heights from a barscan calibration session.

    **Mantid Algorithms used:**
    :ref:`Load <algm-Load-v1>`,

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    barscan_files: str, list
        Path(s) to barscan run file(s). If only one file, it should contain multiple positions of the bar.
        If a list of files, then each file contains the pixel_intensities for every pixel for a fixed position of the
        bar.
    component: str
        Name of the detector panel scanned with the bar. Usually, 'detector1`.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.
    formula: str
        Formula to obtain the position of the bar (Y-coordinate) in the frame of reference located at the sample.
    order: int
        Highest degree for the polynomial that will fit the observed positions of the bar.

    Returns
    -------
    dict
        Dictionary with the following entries:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number associated to the calibration.
        - unit: str, the units for the positions and heights. Set to 'mm' for mili-meters.
        - positions, list, List of Y-coordinate for each pixel.
        - heights, list, List of pixel heights.
    """
    instrument_name, run_number, number_pixels_in_tube, number_tubes = None, None, None, None
    bar_positions = []  # Y-coordinates of the bar for each scan
    # 2D array defining the position of the bar on the detector, in pixel coordinates
    # The first index corresponds to the Y-axis (along each tube), the second to the X-axis (across tubes)
    # Thus, bottom_shadow_pixels[:, 0] indicates bottom shadow pixel coordinates along the very first tubes
    bottom_shadow_pixels = []
    for bar_position, barscan_workspace in barscan_workspace_generator(barscan_files,
                                                                       bar_position_log=bar_position_log):
        if instrument_name is None:
            instrument_name = instrument_standard_name(barscan_workspace)
            run_number = int(SampleLogs(barscan_workspace).single_value('run_number'))
        # Find out the Y-coordinates of the bar in the reference-of-frame located at the sample
        formula_bar_position_inserted = formula.format(y=bar_position)
        bar_positions.append(float(numexpr.evaluate(formula_bar_position_inserted)))
        bottom_shadow_pixels_per_scan = []  # For the current scan, we have one bottom shadow pixel for each tube
        # A TubeCollection is a list of TubeSpectrum objects, representing a physical tube. Here we obtain the
        # list of tubes for the main double-detector-panel.
        # The view 'decreasing X' sort the tubes by decreasing value of their corresponding X-coordinate. In this view,
        # a double detector panel looks like a single detector panel. When looking at the panel standing at the
        # sample, the leftmost tube has the highest X-coordinate, so the 'decreasing X' view orders the tubes
        # from left to right.
        if number_pixels_in_tube is None:
            # We create a tube collection to figure out the pixel indexes for each tube
            collection = TubeCollection(barscan_workspace, component).sorted(view='decreasing X')
            # pixel_indexes is a list of length equal the number of tubes. Each list element is a list containing
            # the pixel indexes for a particular tube.
            pixel_indexes = [tube.spectrum_info_index for tube in collection]
            number_tubes, number_pixels_in_tube = len(collection), len(collection[0])
        pixel_intensities = np.sum(mtd[barscan_workspace].extractY(), axis=1)  # integrated intensity on each pixel
        for pixel_indexes_in_tube in pixel_indexes:  # iterate over each tube, retrieving its pixel indexes
            try:
                # Find the bottom shadow pixel for the current tube and current barscan run
                pixel_intensities_in_tube = pixel_intensities[pixel_indexes_in_tube]
                bottom_shadow_pixels_per_scan.append(find_edges(pixel_intensities_in_tube).bottom_shadow_pixel)
            except Exception:
                # this tube may be malfunctioning for the current barscan
                bottom_shadow_pixels_per_scan.append(INCORRECT_PIXEL_ASSIGNMENT)
        bottom_shadow_pixels.append(bottom_shadow_pixels_per_scan)
    bottom_shadow_pixels = np.array(bottom_shadow_pixels)
    bar_positions = np.array(bar_positions)

    # Deal with corner cases not resolved with the find_edges algorithm
    resolve_incorrect_pixel_assignments(bottom_shadow_pixels, bar_positions)

    if len(bottom_shadow_pixels) <= order:
        raise ValueError(f"There are not enough bar positions to fo a fit with a polynomyal of order {order}.")

    # fit pixel positions for each tube and output in a dictionary
    positions = []
    heights = []
    for tube_index in range(number_tubes):  # iterate over the tubes in the collection
        # Fit the pixel numbers and Y-coordinates of the bar for the current tube with a polynomial
        fit_results = fit_positions(bottom_shadow_pixels[:, tube_index], bar_positions, order=order,
                                    tube_pixels=number_pixels_in_tube)
        # Store the fitted Y-coordinates and heights of each pixel in the current tube
        # Store as lists so that they can be easily serializable
        positions.append(list(fit_results.calculated_positions))  # store with units of mili-meters
        heights.append(list(fit_results.calculated_heights))  # store with units of mili-meters

    return dict(instrument=instrument_name,
                component=component,
                run=run_number,
                unit='mm',
                positions=positions,
                heights=heights)


def resolve_incorrect_pixel_assignments(bottom_shadow_pixels, bar_positions):
    r"""
    # Corner case 1:
    # When the bar approaches the bottom of the tubes, the bottom edge of its shadow blends with the non-functioning
    # pixels at the bottom of the tubes. In this scenario, the bottom edge of the bar is assigned as the first
    # non-functioning pixel at the top of the tube. Thus, we must find out for each tube a sudden jump in the
    # identified bottom shadow pixel.
    # Example: the bottom shadow pixels for the first tube as we change the position of the bar have been identified:
    #    249, 230, 209, 187, 165, 145, 129, 103, 82, 67, 52, 39, 21, 249, 249
    # In the example, the last two bottom pixels are erroneous. They must be set to -1
    #
    # Corner case 2:
    # When the identified bottom pixel in a tube is far from the boundary between the illuminated pixels and the non
    # illuminated pixels, it leaves a few non-illuminated pixels as active pixels. The bottom shadow pixel can then
    # be incorrectly assigned as one of these non-illuminated pixels. The results are outliers pixel indexes.
    # Example: 249 248 246 4 244 242.... Here "4" is in incorrectly assigned pixel index as bottom of the bar
    # We must find these outliers and assign them as incorrect. We use a linear fit to find out outliers
    Parameters
    ----------
    bottom_shadow_pixels

    Returns
    -------

    """
    number_bar_positions, number_tubes = bottom_shadow_pixels.shape
    for tube_index in range(number_tubes):
        # Correct corner case 1
        jump_index = None
        y = bottom_shadow_pixels[:, tube_index]
        for i in range(number_bar_positions - 2, 0, -1):
            if abs(y[i] - y[i-1]) > 10 * max(1, abs(y[i+1] - y[i])):  # value/factor of 10 selected as threshold
                jump_index = i
                break
        if jump_index is not None:
            y[jump_index:] = INCORRECT_PIXEL_ASSIGNMENT
        # Correct corner case 2
        y = bottom_shadow_pixels[:, tube_index]
        x = bar_positions[y != INCORRECT_PIXEL_ASSIGNMENT]
        coefficients = np.polynomial.polynomial.polyfit(x, y[y != INCORRECT_PIXEL_ASSIGNMENT], 1)
        y_fitted = np.polynomial.polynomial.polyval(bar_positions, coefficients)
        residuals = np.abs(y - y_fitted)
        y[residuals > np.average(residuals) + 2 * np.std(residuals)] = INCORRECT_PIXEL_ASSIGNMENT


# Note: a CPP Mantid algorithm similar to ApplyCalibration would be much faster
def apply_barscan_calibration(input_workspace, calibration, output_workspace=None):
    r"""
    Update the pixel positions (Y-coordinate only) and pixel heights of a double-panel in an input workspace.

    **Mantid Algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        Input workspace containing the original pixel positions and heights
    calibration: dict
        Dictionary with the following entries:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number associated to the calibration.
        - unit: str, the units for the positions and heights. Set to 'mm' for mili-meters.
        - positions, list, List of Y-coordinate for each pixel.
        - heights, list, List of pixel heights.
    output_workspace: str
        Name of the workspace containing the updated pixel positions and pixel heights. If :py:obj:`None`, the name of
        ``input_workspace`` is used, therefore modifiying the input workspace. If not :py:obj:`None`, then a clone
         of ``input_workspace`` is produced, but with updated pixel positions and heights.
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)  # pixel positions and heights to be updated for the input workspace
    else:
        # A new workspace identical to the input workspace except in regards to pixel  positions and heights.
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)

    # 2D array of pixel Y-coordinates. The first index is the tube-index
    pixel_positions = calibration['positions']
    factor = 1.e-03 if calibration['unit'] == 'mm' else 1.0  # from mili-meters to meters
    # A TubeCollection is a list of TubeSpectrum objects, representing a physical tube. Here we obtain the
    # list of tubes for the main double-detector-panel.
    # The view 'decreasing X' sort the tubes by decreasing value of their corresponding X-coordinate. In this view,
    # a double detector panel looks like a single detector panel. When looking at the panel standing at the
    # sample, the leftmost tube has the highest X-coordinate, so the 'decreasing X' view orders the tubes
    # from left to right.
    collection = TubeCollection(output_workspace, calibration['component']).sorted(view='decreasing X')
    for tube_index, tube in enumerate(collection):
        if True in np.isnan(pixel_positions[tube_index]):  # this tube was not calibrated
            continue
        # update Y-coord of pixels in this tube. "factor" ensures the units are in meters
        tube.pixel_y = [factor * y for y in pixel_positions[tube_index]]
        tube.pixel_heights = [factor * h for h in calibration['heights'][tube_index]]


def calculate_apparent_tube_width(flood_input, component='detector1', load_barscan_calibration=True):
    r"""
    Determine the tube width most efficient for detecting neutrons. An effective tube (or pixel) diameter is
    determined for tubes in the front panel, and likewise for the tubes in the back panel.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    flood_input: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Path to flood run, flood workspace name, or flood workspace object.
    component: str
        Name of the instrument component containing the detector array consisting of two parallel panels of tubes.
    load_barscan_calibration: bool
        Load pixel positions and heights from the pixel-calibrations database appropriate to ``input_workspace``. If
        py:obj:`False`, then the pixel positions and heigths will be those of ``input_workspace``.

    **Mantid algorithms used:**
        :ref:`DeleteWorkspaces <algm-DeleteWorkspaces-v1>`,
        :ref:`Integration <algm-Integration-v1>`,
        :ref:`MaskDetectors <algm-MaskDetectors-v1>`,
        :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
        :ref:`ReplaceSpecialValues <algm-ReplaceSpecialValues-v1>`,

    Returns
    -------
    Returns
    -------
    dict
        Dictionary containing the following keys:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number of ``input_workspace``.
        - unit: str, the units for the tube widths. Set to 'mm' for mili-meters.
        - widths, list, A two-item list containing the apparent widths for the front and back tubes.
    """
    # Determine the type of input for the flood data
    if file_exists(flood_input):
        input_workspace = unique_workspace_dundername()
        Load(Filename=flood_input, OutputWorkspace=input_workspace)
    else:
        input_workspace = flood_input

    integrated_intensities = unique_workspace_dundername()
    Integration(InputWorkspace=input_workspace, OutputWorkspace=integrated_intensities)

    # Mask non-finite pixel_intensities (nan, inf). They can't be used in the calculation.
    #
    # Replace non-finite pixel_intensities with a value of -1
    ReplaceSpecialValues(InputWorkspace=integrated_intensities, OutputWorkspace=integrated_intensities,
                         NanValue=-1, NanError=-1, InfinityValue=-1, InfinityError=-1)
    # Mask detectors with negative pixel_intensities
    mask_workspace = unique_workspace_dundername()
    MaskDetectorsIf(InputWorkspace=integrated_intensities, Operator='Less', Value=0., OutputWorkspace=mask_workspace)
    MaskDetectors(Workspace=integrated_intensities, MaskedWorkspace=mask_workspace)

    # Update pixel positions and heights with the appropriate calibration, if so requested.
    if load_barscan_calibration is True:
        run_number = SampleLogs(input_workspace).single_value('run_number')
        calibration = load_calibration(instrument_standard_name(input_workspace), run=run_number, component=component)
        apply_barscan_calibration(input_workspace, calibration)

    # Calculate the count density for each tube. Notice that if the whole tube is masked, then the associated
    # intensity is stored as nan.
    #
    # Sort the tubes according to the X-coordinate in decreasing value. This is the order when sitting on the
    # sample and iterating over the tubes "from left to right"
    collection = TubeCollection(integrated_intensities, 'detector1').sorted(view='decreasing X')
    count_densities = list()
    for tube in collection:
        weighted_intensities = tube.readY.ravel() / tube.pixel_heights
        d = np.mean(weighted_intensities[~tube.isMasked])
        count_densities.append(d)
    count_densities = np.array(count_densities)  # is convenient to cast densities into a numpy array data structure.

    # Determine the count densities per panel and for the whole detector array.
    # We must be careful to pick only tubes with finite densities (avoid 'nan')
    average_count_density = np.mean(count_densities[np.isfinite(count_densities)])
    front_count_density = np.mean(count_densities[::2][np.isfinite(count_densities[::2])])  # front tubes, even indexes
    back_count_density = np.mean(count_densities[1::2][np.isfinite(count_densities[1::2])])  # back tubes, odd indexes

    # Determine the front and back pixel widths
    nominal_width = collection[0][0].width  # width of the first pixel in the first tube
    front_width = (front_count_density / average_count_density) * nominal_width
    back_width = (back_count_density / average_count_density) * nominal_width

    DeleteWorkspaces(integrated_intensities, mask_workspace)

    return dict(instrument=instrument_standard_name(input_workspace),
                component=component,
                run=SampleLogs(input_workspace).single_value('run_number'),
                unit='mm',
                widths=(1000 * front_width, 1000 * back_width))


# Note: a CPP Mantid algorithm similar to ApplyCalibration would be much faster
def apply_apparent_tube_width(input_workspace, calibration, output_workspace=None):
    r"""
    Update the pixel widths with effective tube widths.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Input workspace, usually a flood run.
    calibration: dict
        Dictionary with the following required entries:
        - instrument, str, Name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - unit: str, the units for the positions and heights. Usually 'mm' for mili-meters.
        - widths, list, A two-item list containing the apparent widths for the front and back tubes.
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of ``input_workspace`` is used, thus
        calibrating the pixel widths of the input workspace.

    **Mantid algorithms used:**
        :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,

    Returns
    -------
    ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    else:
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)

    collection = TubeCollection(output_workspace, calibration['component']).sorted(view='decreasing X')
    factor = 1.e-03 if calibration['unit'] == 'mm' else 1.0
    front_width,  back_width = factor * np.array(calibration['widths'])

    for tube_index, tube in enumerate(collection):
        tube.pixel_widths = front_width if tube_index % 2 == 0 else back_width  # front tubes have an even index

    return mtd[output_workspace]


"""
Below are functions that combine barscan and flat-field calibrations. These are functions exposed to the public
interfaces.
"""


def calculate_pixel_calibration(barscan_files, flood_file, component='detector1', sample_log='dcal',
                                formula='565 - {dcal}', save_to_database=False):
    r"""
    Calculate pixel positions (only Y-coordinae), pixel heights, and tube widths most efficient for detecting
    neutrons.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    barscan_files: list
        Paths to barscan run files. Each file is a nexus file containing the pixel_intensities for every pixel for a
        given position of the bar.
    flood_file: str
        Path of a flat-field or flood file
    component: str
        Name of the detector panel scanned with the bar. Usually, 'detector1`.
    sample_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.
    formula: str
        Formula to obtain the position of the bar (Y-coordinate) in the frame of reference located at the sample.
    save_to_database: bool
        Option to save the result to the database of pixel-calibrations. Calibrations can later be retrieved and
        applied to workspaces in order to update their pixel positions, heights, and widths.

    Returns
    -------
    dict
        Dictionary with the following entries:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number associated to the calibration.
        - unit: str, the units for the positions and heights. Set to 'mm' for mili-meters.
        - positions, list, List of Y-coordinate for each pixel.
        - heights, list, List of pixel heights.
        - widths, list, A two-item list containing the apparent widths for the front and back tubes.
    """
    # First find out pixel positions and heights using the bar scan files
    barscan_calibration = calculate_barscan_calibration(barscan_files, component=component,
                                                        bar_position_log=sample_log, formula=formula, order=5)
    # Next find out the tube widths using the calibration for pixel positions and heights
    flood_workspace = unique_workspace_dundername()
    Load(flood_file, OutputWorkspace=flood_workspace)
    apply_barscan_calibration(flood_workspace, barscan_calibration)  # update pixel positions and heights
    flood_calibration = calculate_apparent_tube_width(flood_workspace, component='detector1',
                                                      load_barscan_calibration=False)
    # Merge both calibrations
    calibration = barscan_calibration
    calibration['run'] = max(barscan_calibration['run'], flood_calibration['run'])
    calibration['widths'] = flood_calibration['widths']

    if save_to_database is True:
        save_calibration(calibration)

    return calibration


def load_and_apply_pixel_calibration(input_workspace, output_workspace=None, components=['detector1'],
                                     database=None):
    r"""
    Update pixel positions and heights, as well as front and back tube widths, for a specific
    double panel detector array.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    **Mantid algorithms used:**
        :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Input workspace to update
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of ``input_workspace`` is used, thus
        calibrating the pixels of the input workspace.
    components: list
    database: str
        Path to database. If :py:obj:`None`, the default database is used.
    """
    if output_workspace is None:
        output_workspace = input_workspace
    else:
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)

    instrument = instrument_enum_name(output_workspace)
    run_number = SampleLogs(output_workspace).run_number.value
    for component in components:
        calibration = load_calibration(instrument, run=run_number, component=component, database=database)
        if calibration.get('heights', None) is not None:
            apply_barscan_calibration(output_workspace, calibration)
        if calibration.get('widths', None) is not None:
            apply_apparent_tube_width(output_workspace, calibration)


def visualize_barscan_calibration(input_workspace=None, component='detector1', calibration=None,
                                  output_base_workspace=None):
    r"""
    Creates two workspaces with pixel positions and heights as respective intensities

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        Retrieve the pixel positions and heights from this workspace
    component: str
        If ``input_workspace`` is passed, then retrieve pixel positions and heights for this instrument component.
    calibration: dict
        The output of ~drtsans.pixel_calibration.calculate_barscan_calibration
    output_base_workspace: str
        Basename for the two output workspace. Suffixes _positions and _heigths will be appended. If :py:obj:`None`,
        then name of ``input_workspace`` is used if provided.
    """
    if input_workspace is None and calibration is None:
        raise RuntimeError('Provide either an input workspace or a calibration')

    if output_base_workspace is None:
        if input_workspace is None:
            raise RuntimeError('Please provide and output basename for the output workspaces')
        else:
            output_base_workspace = str(input_workspace)

    # Generate a calibration from the provided input workspace and component
    if input_workspace is not None:
        calibration['instrument'] = instrument_standard_name(input_workspace)
        calibration['unit'] = 'mm'
        calibration['positions'] = list()
        calibration['heights'] = []
        collection = TubeCollection(input_workspace, component=component).sorted(view='decreasing X')
        for tube in collection:
            calibration['positions'].append(list(1000 * tube.pixel_y))
            calibration['heights'].append(list(1000 * tube.pixel_heights))

    #
    for feature in ('positions', 'heights'):
        y = np.array(calibration[feature]).ravel()
        if feature == 'positions':
            y -= min(y)  # minimum position set to zero
        CreateWorkspace(DataX=[0, 1], DataY=y, NSpec=len(y),
                        OutputWorkspace=output_base_workspace + '_' + str(feature))
        LoadInstrument(Workspace=output_base_workspace + '_' + str(feature),
                       InstrumentName=calibration['instrument'], RewriteSpectraMap=True)
