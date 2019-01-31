from __future__ import absolute_import, division, print_function

import numpy as np

from mantid.kernel import logger
from mantid.simpleapi import CloneWorkspace, MaskDetectors


def _interpolate_tube(x, y, e, detectors_masked, polynomial_degree):

    xx = x[~detectors_masked]
    yy = y[~detectors_masked]
    ee = e[~detectors_masked]

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
    # Try this: There is a method in `MatrixWorkspace`, 
    # `getDetectorIDToWorkspaceIndexMap` or something like that.
    # Not exposed to python :(
    spectrum_info = flood_ws.spectrumInfo()
    ws_index = 0
    while spectrum_info.isMonitor(ws_index):
        logger.notice("No monitor: ws_index:  {}".format(ws_index))
        ws_index += 1
    
    first_detector_index = ws_index
    logger.notice("First detector index: ws_index:  {}".format(first_detector_index))

    # Put the data in numpy arrays
    data_y = flood_ws.extractY()
    data_e = flood_ws.extractE()

    # numpy arrays with all the detector IDs
    detector_info = flood_ws.detectorInfo()
    detector_ids = detector_info.detectorIDs()

    number_of_tubes = component.nelements()
    # Lets get the output workspace
    __output_ws = CloneWorkspace(flood_ws)

    for tube_idx in range(number_of_tubes):
        logger.notice("Tube index beeing analised = {}. WS index={}.".format(
            tube_idx, ws_index))

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

        tube_detector_ids = detector_ids[tube_idx*num_detectors_in_tube: num_detectors_in_tube + tube_idx*num_detectors_in_tube]
        print("********** tube_detector_ids **************\n", tube_detector_ids)

            # ws_index-first_detector_index:ws_index-first_detector_index+num_detectors_in_tube]
        # check how many detectors are masked in the tube
        # same size as the tube, True where masked
        detectors_masked = np.array(
            [detector_info.isMasked(int(id)) for id in tube_detector_ids])
        # Count the detectors masked (True)
        num_of_detectors_masked = np.count_nonzero(detectors_masked)
        if num_of_detectors_masked > 0 and \
                (len(tube_detector_ids) - num_of_detectors_masked) > \
                    min_detectors_per_tube:  # number of detectors with values > min_detectors_per_tube
            
            # Let's fit
            y = data_y[ws_index:ws_index+num_detectors_in_tube].flatten()
            e = data_e[ws_index:ws_index+num_detectors_in_tube].flatten()

            y_new, e_new = _interpolate_tube(x, y, e, detectors_masked,
                                             polynomial_degree)

            for detector_id, y_new_value, e_new_value in zip(
                    tube_detector_ids[detectors_masked], y_new[detectors_masked],
                    e_new[detectors_masked]):
                # int because we are iterating numpy array and mantid set needs
                # type int
                ws_index_to_edit = int(detector_id - first_detector_index)
                # logger.notice("{} setY({},{})  setE({},{})".format(
                #     __output_ws.name(), ws_index_to_edit, y_new_value,
                #     ws_index_to_edit, e_new_value))

                __output_ws.setY(ws_index_to_edit, np.array([y_new_value]))
                __output_ws.setE(ws_index_to_edit, np.array([e_new_value]))
                detector_info.setMasked(int(detector_id), True)

        elif num_of_detectors_masked > min_detectors_per_tube:
            logger.error("Skipping tube {}. Too many masked detectors.".format(
                tube_idx))

        ws_index += num_detectors_in_tube

    return __output_ws
