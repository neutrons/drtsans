import numpy as np
import numexpr
import tinydb
import collections

r""" Hyperlinks to mantid algorithms
CloneWorkspace <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
DeleteWorkspaces <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspaces-v1.html>
Integration <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
MaskDetectors <https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html>
MaskDetectorsIf <https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html>
ReplaceSpecialValues <https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html>
"""
from mantid.simpleapi import (CloneWorkspace, DeleteWorkspaces, Integration, Load, MaskDetectors, MaskDetectorsIf,
                              ReplaceSpecialValues)
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


def _consecutive_true_values(values, how_many, reverse=False,
                             message='Could not find pixel index'):
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
    message: str
        Exception message

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
        raise IndexError(message)


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
        pixel intensities along the tube.
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
    # calculate minimum intensity thresholds for tube ends and shaddows
    av = np.average(intensities)
    end_threshold = tube_threshold * av
    shadow_threshold = shadow_threshold * av

    # Find edges of the tube: want at least tube_edge_min_width pixels
    # (starting from the top or bottom of a tube) that have intensities greater
    # than the threshold
    illuminated = [bool(i > end_threshold) for i in intensities]
    bottom_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                            message='Could not find bottom tube edge')
    top_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                         message='Could not find top tube edge',
                                         reverse=True)

    # Find the shadow region: similar to tube edges, but in this case
    # we want shadow_edge_min_width intensities less than the shadow threshold,
    # followed by at least one intensity greater than the threshold
    active = intensities[bottom_pixel: top_pixel+1]
    shadowed = [bool(i < shadow_threshold) for i in active]
    bottom_shadow_pixel = bottom_pixel +\
        _consecutive_true_values(shadowed, shadow_edge_min_width,
                                 message='Could not find bottom shadow edge')

    active = intensities[bottom_shadow_pixel + shadow_edge_min_width:
                         top_pixel + 1]
    illuminated = [bool(i > shadow_threshold) for i in active]
    above_shadow_pixel = bottom_shadow_pixel + shadow_edge_min_width +\
        _consecutive_true_values(illuminated, 1,
                                 message='could not find above shadow region')

    # Check for a faulty tube: we want a certain number of pixels not in the bar shaddow
    active_tube_length = top_pixel - bottom_pixel + 1
    shadow_length = above_shadow_pixel - bottom_shadow_pixel
    if active_tube_length < min_illuminated_length + shadow_length:
        raise RuntimeError('Faulty tube found')

    return dict(bottom_pixel=bottom_pixel, top_pixel=top_pixel,
                bottom_shadow_pixel=bottom_shadow_pixel,
                above_shadow_pixel=above_shadow_pixel)


@namedtuplefy
def fit_positions(edge_pixels, bar_positions, tube_pixels=256, order=5):
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

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - calculated_positions: calculated positions of the pixels
        - calculated_heights: calculated pixel heights
    """
    message_len = 'The positions of the bar and edge pixels have to be the same length'
    assert len(edge_pixels) == len(bar_positions), message_len

    try:
        # fit the bar positions to a 5th degree polynomial in edge_pixels
        coefficients = np.polynomial.polynomial.polyfit(edge_pixels, bar_positions, int(order))
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


def load_calibration(instrument, run=None, component='detector1'):
    r"""
    Load pixel calibration from the database.

    Parameters
    ----------
    instrument:str
        Name of the instrument
    run: int
        Run number to resolve which calibration to use. If py:obj:`None`, then the lates calibration will be used.
    component: str
        Name of the double panel detector array for which the calibration was performed

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
    database = tinydb.TinyDB(database_file[enum_instrument])
    calibrations = database.Query()

    def find_matching_calibration(calibration_type_query):
        matches = database.search(calibration_type_query &
                                  (calibrations.instrument == enum_instrument) &
                                  (calibrations.component == component))
        # Sort by decreasing run number and find appropriate calibration
        matches = sorted(matches, key=lambda d: d[run], reverse=True)
        if run is None:
            return matches[0][1]  # return latest calibration
        else:
            for i, match in enumerate(matches):
                if run > match['run']:
                    return match[i-1]  # return first calibration with a smaller run number than the query run number

    # Find matching barscan calibration (pixel positions and heights)
    calibration = find_matching_calibration(calibrations.heights)
    # Find matching flat-field calibration (pixel widths) and update the calibration dictionary
    calibration['widths'] = find_matching_calibration(calibrations.widths)['widths']

    return calibration


def save_calibration(calibration):
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
    tinydb.TinyDB(database_file[enum_instrument]).insert(dict(calibration))


def calculate_barscan_calibration(barscan_files, component='detector1', sample_log='dcal',
                                  formula='565 - {dcal}', order=5):
    r"""
    Calculate pixel positions (only Y-coordinae) as well as pixel heights from a barscan calibration session.

    **Mantid Algorithms used:**
    :ref:`Load <algm-Load-v1>`,

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    barscan_files: list
        Paths to barscan run files. Each file is a nexus file containing the intensities for every pixel for a given
        position of the bar.
    component: str
        Name of the detector panel scanned with the bar. Usually, 'detector1`.
    sample_log: str
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
    if len(barscan_files) <= order:
        raise ValueError(f"There are not enough files to fo a fit with a polynomyal of order {order}.")
    bar_positions = []  # Y-coordinates of the bar for each scan
    # 2D array defining the position of the bar on the detector, in pixel coordinates
    # The first index corresponds to the Y-axis (along each tube), the second to the X-axis (across tubes)
    # Thus, bottom_shadow_pixels[:, 0] indicates bottom shadow pixel coordinates along the very first tubes
    bottom_shadow_pixels = []
    workspace_name = unique_workspace_dundername()
    # loop over each barscan run, stored in each of the barscan_files
    for filename in barscan_files:
        Load(filename, OutputWorkspace=workspace_name)
        # Find out the Y-coordinates of the bar in the reference-of-frame located at the sample
        formula_dcal_inserted = formula.format(dcal=SampleLogs(workspace_name).find_log_with_units(sample_log, 'mm'))
        bar_positions.append(float(numexpr.evaluate(formula_dcal_inserted)))
        bottom_shadow_pixels_per_scan = []  # For the current scan, we have one bottom shadow pixel for each tube
        # A TubeCollection is a list of TubeSpectrum objects, representing a physical tube. Here we obtain the
        # list of tubes for the main double-detector-panel.
        # The view 'decreasing X' sort the tubes by decreasing value of their corresponding X-coordinate. In this view,
        # a double detector panel looks like a single detector panel. When looking at the panel standing at the
        # sample, the leftmost tube has the highest X-coordinate, so the 'decreasing X' view orders the tubes
        # from left to right.
        collection = TubeCollection(workspace_name, component).sorted(view='decreasing X')
        for tube in collection:  # iterate over each tube in the collection of tubes
            try:
                # Find the bottom shadow pixel for the current tube and current barscan run
                bottom_shadow_pixels_per_scan.append(find_edges(tube.readY.ravel()).bottom_shadow_pixel)
            except Exception:
                bottom_shadow_pixels_per_scan.append(np.nan)  # this tube may be malfunctioning for the current barscan
        bottom_shadow_pixels.append(bottom_shadow_pixels_per_scan)
    bottom_shadow_pixels = np.array(bottom_shadow_pixels)

    # fit pixel positions for each tube and output in a dictionary
    positions = []
    heights = []
    collection = TubeCollection(workspace_name, component)
    number_pixels_in_tube = len(collection[0])  # length of first tube
    for tube_index in range(len(collection)):  # iterate over the tubes in the collection
        # Fit the pixel numbers and Y-coordinates of the bar for the current tube with a polynomial
        fit_results = fit_positions(bottom_shadow_pixels[:, tube_index], bar_positions, order=order,
                                    tube_pixels=number_pixels_in_tube)
        # Store the fitted Y-coordinates and heights of each pixel in the current tube
        positions.append(fit_results.calculated_positions)  # store with units of mili-meters
        heights.append(fit_results.calculated_heights)  # store with units of mili-meters

    return dict(instrument=instrument_standard_name(workspace_name),
                component=component,
                run=SampleLogs(workspace_name).single_value('run_number'),
                unit='mm',
                positions=positions,
                heights=heights)


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

    # Mask non-finite intensities (nan, inf). They can't be used in the calculation.
    #
    # Replace non-finite intensities with a value of -1
    ReplaceSpecialValues(InputWorkspace=integrated_intensities, OutputWorkspace=integrated_intensities,
                         NanValue=-1, NanError=-1, InfinityValue=-1, InfinityError=-1)
    # Mask detectors with negative intensities
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
        Paths to barscan run files. Each file is a nexus file containing the intensities for every pixel for a given
        position of the bar.
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
    barscan_calibration = calculate_barscan_calibration(barscan_files, component=component, sample_log=sample_log,
                                                        formula=formula, order=5)
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


def load_and_apply_pixel_calibration(input_workspace, output_workspace=None, components=['detector1']):
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
        Names of the double panel detector arrays.
    """
    if output_workspace is None:
        output_workspace = input_workspace
    else:
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)

    instrument = instrument_enum_name(output_workspace)
    run_number = SampleLogs(output_workspace).single_value('run_number')
    for component in components:
        calibration = load_calibration(instrument, run=run_number, component=component)
        apply_barscan_calibration(output_workspace, calibration)
        apply_apparent_tube_width(output_workspace, calibration)
