from __future__ import absolute_import, division, print_function

import numpy as np

from mantid.kernel import logger
from mantid.simpleapi import CloneWorkspace, MaskDetectors


def _interpolate_tube(x, y, e, pixels_masked, polynomial_degree):

    xx = x[~pixels_masked]
    yy = y[~pixels_masked]
    ee = e[~pixels_masked]

    polynomial_coeffs_y = np.polyfit(xx, yy, polynomial_degree)
    polynomial_coeffs_e = np.polyfit(xx, ee, polynomial_degree)

    y_new = np.polyval(polynomial_coeffs_y, x)
    e_new = np.polyval(polynomial_coeffs_e, x)

    return y_new, e_new


def interpolate_mask(flood_ws, polynomial_degree=1,
                     component_name='detector1', min_pixels_per_tube=90):
    '''Interpolates over the mask (usually the beamstop and the tube ends)
    Assumptions:
    - Monitors are in the beginning of the workspace
    - IDs of the detectors are contiguous

    Parameters
    ----------
    flood_ws : MatrixWorkspace
        input workspace. Normally a sensitivity file.
    polynomial_degree : int, optional
        Polynomial degree for the interpolation (the default is 1, which is a
        linear fitting). A value of 0 calculates de average.
    component_name : str, optional
        Component name to  (the default is 'detector1', which is the main
        detector)
    min_pixels_per_tube : int, optional
        Minimum pixels existing in the tube to fit (the default is 50)

    TODO: Average when polynomial_degree == 2

    Returns
    -------
    MatrixWorkspace
        The interpolated workspace
    '''

    # Get the main detector: the masks are in the mask workspace!!!
    instrument = flood_ws.getInstrument()
    component = instrument.getComponentByName(component_name)
    if not component:
        logger.error("Component not valid! {}".format(component_name))
        return

    # find out in which WS index starts the monitor
    # Try this: There is a method in `MatrixWorkspace`, `getDetectorIDToWorkspaceIndexMap` or something like that.
    spectrumInfo = flood_ws.spectrumInfo()
    ws_index = 0
    while spectrumInfo.isMonitor(ws_index):
        ws_index += 1
    first_detector_index = ws_index
    # Put the data in numpy arrays
    data_y = flood_ws.extractY()
    data_e = flood_ws.extractE()

    # numpy arrays with all the detector IDs
    detectorInfo = flood_ws.detectorInfo()
    detectorIDs = detectorInfo.detectorIDs()

    #
    spectrumInfo = flood_ws.spectrumInfo()

    number_of_tubes = component.nelements()

    __interpolated_ws = CloneWorkspace(flood_ws)

    # Iterate through tubes: from 0 to number_of_tubes
    num_pixels_in_tube = 0
    x = []  # this is going to be the tube Y coordinates
    for tube_idx in range(number_of_tubes):
        logger.notice("Tube index beeing analised: {}".format(tube_idx))

        if component[0].nelements() <= 1:
            # Handles EQSANS
            tube = component[tube_idx][0]
        else:
            # Handles Biosans/GPSANS
            tube = component[tube_idx]

        # first tube operations
        if tube_idx == 0:
            num_pixels_in_tube = tube.nelements()
            # x is the Y coordinates of every pixel in the tube
            x = np.array([tube[i].getPos()[1] for i in range(
                num_pixels_in_tube)])

        # check how many pixels are masked in the tube
        tubeIDs = detectorIDs[ws_index:ws_index+num_pixels_in_tube]
        # same size as the tube, True where masked
        pixels_masked = np.array(
            [detectorInfo.isMasked(int(id)) for id in tubeIDs])

        num_of_pixels_masked = np.count_nonzero(pixels_masked)
        if num_of_pixels_masked > 0 and \
                num_of_pixels_masked <= min_pixels_per_tube:

            # Let's fit
            y = data_y[ws_index:ws_index+num_pixels_in_tube].flatten()
            e = data_e[ws_index:ws_index+num_pixels_in_tube].flatten()

            y_new, e_new = _interpolate_tube(x, y, e, pixels_masked,
                                             polynomial_degree)

            for detectorID, y_new_value, e_new_value in zip(
                    tubeIDs[~pixels_masked], y_new[~pixels_masked],
                    e_new[~pixels_masked]):
                
                detector_id_to_replace = detectorID - first_detector_index
                
                __interpolated_ws.setY(int(detector_id_to_replace),
                                       np.array([y_new_value]))
                __interpolated_ws.setE(int(detector_id_to_replace),
                                       np.array([e_new_value]))

        elif num_of_pixels_masked > min_pixels_per_tube:
            logger.error("Skipping tube {}. Too many masked pixels.".format(
                tube_idx))

        ws_index += num_pixels_in_tube

    return __interpolated_ws
