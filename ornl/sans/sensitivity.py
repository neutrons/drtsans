from __future__ import absolute_import, division, print_function

import numpy as np

from mantid.kernel import logger
from mantid.simpleapi import CloneWorkspace, MaskDetectors


def _interpolate_tube(x, y, e, pixelsMasked, polynomial_degree):

    xx = x[~pixelsMasked]
    yy = y[~pixelsMasked]
    ee = e[~pixelsMasked]

    polynomial_coeffs_y = np.polyfit(xx, yy, polynomial_degree)
    polynomial_coeffs_e = np.polyfit(xx, ee, polynomial_degree)

    y_new = np.polyval(polynomial_coeffs_y, x)
    e_new = np.polyval(polynomial_coeffs_e, x)

    return y_new, e_new


def patch_mask(data_ws, mask_ws, polynomial_degree=1,
               component_name='detector1', min_pixels_per_tube=90):
    '''Interpolates over the mask (usually the beamstop and the tube ends)
    Assumptions:
    - Monitors are in the beginning of the workspace
    - IDs of the detectors are contiguous

    Parameters
    ----------
    data_ws : MatrixWorkspace
        input workspace
    mask_ws : ?
        The masked WS where, usually, the beamstop and the tube ends are.
    polynomial_degree : int, optional
        Polynomial degree for the interpolation (the default is 1, which is a
        linear fitting). A value of 0 calculates de average.
    component_name : str, optional
        Component name to  (the default is 'detector1', which is the main
        detector)
    min_pixels_per_tube : int, optional
        Minimum pixels existing in the tube to fit (the default is 50)


    '''

    # In case there are masks (e.g. bad tubes) in the data, copy them to
    # the mask.
    # Do we need this???
    # Usually the mask is drawn and saved in the instrument view
    MaskDetectors(Workspace=mask_ws, MaskedWorkspace=data_ws)

    # Get the main detector: the masks are in the mask workspace!!!
    instrument = mask_ws.getInstrument()
    component = instrument.getComponentByName(component_name)
    if not component:
        logger.error("Component not valid! {}".format(component_name))
        return

    # find out in which WS index starts the monitor
    spectrumInfo = data_ws.spectrumInfo()
    ws_index = 0
    while spectrumInfo.isMonitor(ws_index):
        ws_index += 1
    first_detector_index = ws_index
    # Put the data in numpy arrays
    data_y = data_ws.extractY()
    data_e = data_ws.extractE()

    # numpy arrays with all the detector IDs
    detectorInfo = mask_ws.detectorInfo()
    detectorIDs = detectorInfo.detectorIDs()

    #
    spectrumInfo = data_ws.spectrumInfo()

    number_of_tubes = component.nelements()

    output_ws = CloneWorkspace(data_ws)

    # Iterate through tubes: from 0 to number_of_tubes
    num_pixels_in_tube = 0
    x = []  # this is going to be the tube Y coordinates
    for tube_idx in range(number_of_tubes):
        if component[0].nelements() <= 1:
            # Handles EQSANS
            tube = component[tube_idx][0]
        else:
            # Handles Biosans/GPSANS
            tube = component[tube_idx]

        # first tube operations
        if tube_idx == 0:
            num_pixels_in_tube = tube.nelements()
            # x is the Y coordinates of every pixel inn the tube
            x = np.array([tube[i].getPos()[1] for i in range(
                num_pixels_in_tube)])

        # check how many pixels are masked in the tube
        tubeIDs = detectorIDs[ws_index:ws_index+num_pixels_in_tube]
        # same size as the tube, True where masked
        pixelsMasked = np.array(
            [detectorInfo.isMasked(int(id)) for id in tubeIDs])

        num_of_pixels_masked = np.count_nonzero(pixelsMasked)
        if num_of_pixels_masked > 0 and \
                num_of_pixels_masked <= min_pixels_per_tube:

            # Let's fit
            y = data_y[ws_index:ws_index+num_pixels_in_tube].flatten()
            e = data_e[ws_index:ws_index+num_pixels_in_tube].flatten()

            y_new, e_new = _interpolate_tube(x, y, e, pixelsMasked,
                                             polynomial_degree)

            for detectorID, y_new_value, e_new_value in zip(
                    tubeIDs[~pixelsMasked], y_new[~pixelsMasked],
                    y_new[~pixelsMasked]):

                output_ws.setY(detectorID-first_detector_index, [y_new_value])
                output_ws.setE(detectorID-first_detector_index, [e_new_value])

        elif num_of_pixels_masked > min_pixels_per_tube:
            logger.error("Skipping tube {}. Too many masked pixels.".format(
                tube_idx))

        ws_index += num_pixels_in_tube

    return output_ws
