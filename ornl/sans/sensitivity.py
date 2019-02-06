from __future__ import absolute_import, division, print_function

import numpy as np

from mantid.kernel import logger, Property
from mantid.simpleapi import CloneWorkspace


def _interpolate_tube(x, y, e, detectors_masked, detectors_inf,
                      polynomial_degree):

    xx = x[~detectors_masked & ~detectors_inf]
    yy = y[~detectors_masked & ~detectors_inf]
    ee = e[~detectors_masked & ~detectors_inf]

    polynomial_coeffs_y = np.polyfit(xx, yy, polynomial_degree)
    polynomial_coeffs_e = np.polyfit(xx, ee, polynomial_degree)

    y_new = np.polyval(polynomial_coeffs_y, x)
    e_new = np.polyval(polynomial_coeffs_e, x)

    return y_new, e_new


def interpolate_mask(flood_ws, polynomial_degree=1,
                     component_name='detector1', min_detectors_per_tube=50):
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
    min_detectors_per_tube : int, optional
        Minimum detectors with a value existing in the tube to fit
        (the default is 50)

    TODO: Average when polynomial_degree == 0

    Returns
    -------
    MatrixWorkspace
        The interpolated workspace
    '''

    assert flood_ws.getNumberBins() == 1, "This only supports integrated WS"

    # Get the main detector: the masks are in the mask workspace!!!
    instrument = flood_ws.getInstrument()
    component = instrument.getComponentByName(component_name)
    if not component:
        logger.error("Component not valid! {}".format(component_name))
        return

    # find out in which WS index starts the monitor
    # Try this: There is a method in `MatrixWorkspace`,
    # `getDetectorIDToWorkspaceIndexMap` or something like that.
    # Not exposed to python :(
    spectrum_info = flood_ws.spectrumInfo()
    ws_index = 0
    while spectrum_info.isMonitor(ws_index):
        logger.notice("No monitor: ws_index:  {}".format(ws_index))
        ws_index += 1

    first_detector_index = ws_index

    # Put the data in numpy arrays
    data_y = flood_ws.extractY()
    data_e = flood_ws.extractE()

    # numpy arrays with all the detector IDs
    detector_info = flood_ws.detectorInfo()
    detector_ids = detector_info.detectorIDs()

    detector_id_to_index = {id: index for (id, index) in
                            zip(detector_ids,
                                range(first_detector_index,
                                      flood_ws.getNumberHistograms()))}

    number_of_tubes = component.nelements()
    # Lets get the output workspace
    __output_ws = CloneWorkspace(flood_ws)
    detector_info_output_ws = __output_ws.detectorInfo()

    for tube_idx in range(number_of_tubes):

        if component[0].nelements() <= 1:
            # Handles EQSANS
            tube = component[tube_idx][0]
        else:
            # Handles Biosans/GPSANS
            tube = component[tube_idx]

        # first tube operations
        if tube_idx == 0:
            num_detectors_in_tube = tube.nelements()
            # x is the Y coordinates of every detector in the tube
            x = np.array([tube[i].getPos()[1] for i in range(
                num_detectors_in_tube)])

        tube_detector_ids = detector_ids[tube_idx *
                                         num_detectors_in_tube:
                                         num_detectors_in_tube +
                                         tube_idx*num_detectors_in_tube]

        # check how many detectors are masked in the tube
        # same size as the tube, True where masked
        detectors_masked = np.array(
            [detector_info.isMasked(detector_id_to_index[id])
             for id in tube_detector_ids])
        # Those are the detectors marked "inf" (EMPTY_DBL)
        detectors_inf = np.array([flood_ws.readY(
            detector_id_to_index[id])[0] == Property.EMPTY_DBL
            for id in tube_detector_ids])

        # Count the detectors masked (True)
        num_of_detectors_masked = np.count_nonzero(detectors_masked)
        num_of_detectors_inf = np.count_nonzero(detectors_inf)
        if num_of_detectors_masked > 0 and \
                (len(tube_detector_ids) - num_of_detectors_masked -
                 num_of_detectors_inf) > min_detectors_per_tube:
                # number of detectors with values > min_detectors_per_tube
            # Let's fit
            y = data_y[ws_index:ws_index+num_detectors_in_tube].flatten()
            e = data_e[ws_index:ws_index+num_detectors_in_tube].flatten()

            y_new, e_new = _interpolate_tube(x, y, e, detectors_masked,
                                             detectors_inf, polynomial_degree)
            # Set output workspace with interpolated values
            for detector_id, y_new_value, e_new_value in zip(
                    tube_detector_ids[detectors_masked],
                    y_new[detectors_masked], e_new[detectors_masked]):

                ws_index_to_edit = detector_id_to_index[detector_id]

                __output_ws.setY(ws_index_to_edit, np.array([y_new_value]))
                __output_ws.setE(ws_index_to_edit, np.array([e_new_value]))
                detector_info_output_ws.setMasked(ws_index_to_edit, False)

            # mask pixels in WS Out where detectors_inf
            for detector_id in tube_detector_ids[detectors_inf]:
                ws_index_to_edit = detector_id_to_index[detector_id]
                detector_info_output_ws.setMasked(ws_index_to_edit, True)
                __output_ws.setY(ws_index_to_edit, np.array([1]))
                __output_ws.setE(ws_index_to_edit, np.array([0]))

        elif num_of_detectors_masked > min_detectors_per_tube:
            logger.error("Skipping tube {}. Too many masked detectors.".format(
                tube_idx))

        ws_index += num_detectors_in_tube

    return __output_ws


def inf_value_to_mask(ws):
    """In case the input workspace has infinite values (EMPTY_DBL)
    Mask these values

    Parameters
    ----------
    ws : [type]
        [description]

    """

    assert ws.getNumberBins() == 1, "This only supports integrated WS"

    detector_info = ws.detectorInfo()
    for i in range(ws.getNumberHistograms()):
        if ws.readY(i)[0] == Property.EMPTY_DBL:
            detector_info.setMasked(i, True)
            ws.setY(i, np.array([1]))
            ws.setE(i, np.array([0]))
