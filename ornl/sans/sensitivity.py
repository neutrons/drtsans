from __future__ import absolute_import, division, print_function

import numpy as np

from collections import OrderedDict
from mantid.kernel import Property, logger
from mantid.simpleapi import CloneWorkspace
from ornl.settings import unique_workspace_name


class Detector(object):

    def __init__(self, workspace, component_name):
        self._workspace = workspace
        self.n_tubes = None
        self.n_pixels_per_tube = None
        self.first_det_id = None
        self.last_det_id = None
        self.data_y = None
        self.data_e = None
        self.tube_ws_indices = None
        self._current_start_ws_index = None 
        self._current_stop_ws_index = None

        self._detector_details(component_name)
        self._detector_id_to_ws_index()
        self._extract_data()
        self._set_tube_ws_indices()

    def _detector_details(self, component_name):
        '''Sets the details of the detector:
        n_tubes
        n_pixels_per_tube
        first_det_id
        last_det_id

        Parameters
        ----------
        component_name : string
        
        '''

        i = self._workspace.getInstrument()
        component = i.getComponentByName(component_name)
        self.n_tubes = component.nelements()

        if component[0].nelements() == 1:
            # Handles EQSANS
            self.n_pixels_per_tube = component[0][0].nelements()
            self.first_det_id = component[0][0][0].getID()
            self.last_det_id = component[
                self.n_tubes - 1][0][self.n_pixels_per_tube-1].getID()
        else:
            # Handles BioSANS/GPSANS
            self.n_pixels_per_tube = component[0].nelements()
            self.first_det_id = component[0][0].getID()
            self.last_det_id = component[
                self.n_tubes - 1][self.n_pixels_per_tube-1].getID()

    def _detector_id_to_ws_index(self):
        ''' Converts the detector ID of every pixel in workspace indices.
        Sets `detector_id_to_ws_index`.
        '''

        spectrum_info = self._workspace.spectrumInfo()
        detector_id_to_index = OrderedDict(
            (self._workspace.getSpectrum(
                ws_index).getDetectorIDs()[0], ws_index)
            for ws_index in range(self._workspace.getNumberHistograms())
            if not spectrum_info.isMonitor(ws_index))
        self.detector_id_to_ws_index = detector_id_to_index

    def _extract_data(self):
        '''Sets:
        data_y
        data_e
        '''

        self.data_y = self._workspace.extractY()
        self.data_e = self._workspace.extractE()

    def _set_tube_ws_indices(self):
        '''Seta `tube_ws_indices` as an iterator for the next round of tube
        indices beeing callled
        '''

        tube_ws_indices = []
        for tube_idx in range(self.n_tubes):
            first_det_id = self.first_det_id + tube_idx*self.n_pixels_per_tube
            last_det_id = first_det_id + self.n_pixels_per_tube

            first_ws_index = self.detector_id_to_ws_index[first_det_id]
            last_ws_index = self.detector_id_to_ws_index[last_det_id]

            tube_ws_indices.append((first_ws_index, last_ws_index))
        self._tube_ws_indices = iter(tube_ws_indices)
    
    def next_tube(self):
        '''Calls the next iteration of the iterator `self._tube_ws_indices`
        and sets the member variables:
        _current_start_ws_index 
        _current_stop_ws_index 
        '''

        self._current_start_ws_index, self._current_stop_ws_index = next(
            self._tube_ws_indices)
    
    def get_current_ws_indices(self):
        return (self._current_start_ws_index, self._current_stop_ws_index)

    def ws_data(self):
        '''Returns the current tube data and error
        
        Returns
        -------
        Tuple of arrays
            Y and Error
        '''

        start_ws_index = self._current_start_ws_index
        stop_ws_index = self._current_stop_ws_index
        return (self.data_y[start_ws_index:stop_ws_index].flatten(),
                self.data_e[start_ws_index:stop_ws_index].flatten())

    def pixels_masked(self):
        '''Returns an array of booleans for this tube
        where the detector is masked
        '''
        start_ws_index = self._current_start_ws_index
        stop_ws_index = self._current_stop_ws_index
        detector_info = self._workspace.detectorInfo()
        return np.array([detector_info.isMasked(idx) for idx in range(
            start_ws_index, stop_ws_index)])
    
    def pixels_infinite(self):
        '''Returns an array of booleans for this tube
        where the pixel count is EMPTY_DBL
        '''
        start_ws_index = self._current_start_ws_index
        stop_ws_index = self._current_stop_ws_index

        return np.array([self._workspace.readY(idx)[0] == Property.EMPTY_DBL
                         for idx in range(start_ws_index, stop_ws_index)])
        


def _interpolate_tube(x, y, e, detectors_masked, detectors_inf,
                      polynomial_degree):
    '''Interpolates a tube.
    Calculates the fitting for non masked pixels and non dead pixels

    Parameters
    ----------
    x : np.array
        linear position of the pixel
    y : np.array
        values of pixels
    e : np.array
        Error on the value of a pixel
    detectors_masked : np.array
        Same size as x,y,z, True where pixel is masked
    detectors_inf : np.array
        Same size as x,y,z, True where pixel is dead

    Returns
    -------
    Tuple (y,error)
    '''

    # Get the values to calculate the fit
    xx, yy, ee = [arr[~detectors_masked & ~detectors_inf] for arr in (x, y, e)]

    polynomial_coeffs, cov_matrix = np.polyfit(
        xx, yy, polynomial_degree, w=1/ee, cov=True)
    y_new = np.polyval(polynomial_coeffs, x)
    # Errors in the least squares is the sqrt of the covariance matrix
    # (correlation between the coefficients)
    e_coeffs = np.sqrt(np.diag(cov_matrix))
    # errors of the polynomial
    e_new = np.sqrt([np.add.reduce(
        [(e_coeff*n**i)**2 for i, e_coeff in enumerate(e_coeffs)]) for n in x])

    return y_new, e_new


def interpolate_mask(flood_ws, polynomial_degree=1,
                     component_name='detector1', min_detectors_per_tube=50):
    '''Interpolates over the mask (usually the beamstop and the tube ends)

    Parameters
    ----------
    flood_ws : MatrixWorkspace
        input workspace. Normally a sensitivity file.
    polynomial_degree : int, optional
        Polynomial degree for the interpolation (the default is 1, which is a
        linear fitting).
    component_name : str, optional
        Component name to  (the default is 'detector1', which is the main
        detector)
    min_detectors_per_tube : int, optional
        Minimum detectors with a value existing in the tube to fit. Only fits
        tubes with at least `min_detectors_per_tube` (the default is 50).

    TODO: Refactor this crap?

    Returns
    -------
    MatrixWorkspace
        The interpolated workspace
    '''

    # Sanity check
    assert flood_ws.getNumberBins() == 1, "This only supports integrated WS"

    # Get the main detector: the masks are in the mask workspace!!!
    instrument = flood_ws.getInstrument()
    component = instrument.getComponentByName(component_name)

    detector_id_to_index = _get_detector_id_to_ws_index(flood_ws)
    detector_ids = list(detector_id_to_index.keys())
    ws_indices = list(detector_id_to_index.values())
    ws_index = ws_indices[0]

    # Put the data in numpy arrays
    data_y = flood_ws.extractY()
    data_e = flood_ws.extractE()

    number_of_tubes = component.nelements()
    # Lets get the output workspace
    output_ws = CloneWorkspace(
        flood_ws, OutputWorkspace=unique_workspace_name(
            prefix="__sensitivity_"))

    detector_info_input_ws = flood_ws.detectorInfo()
    detector_info_output_ws = output_ws.detectorInfo()

    for tube_idx in range(number_of_tubes):
        tube = _get_tube(component, tube_idx)

        # first tube operations
        if tube_idx == 0:
            num_detectors_in_tube = tube.nelements()
            # x is the Y coordinates of every detector in the tube
            x = np.array([tube[i].getPos()[1] for i in range(
                num_detectors_in_tube)])

        tube_detector_ids = detector_ids[
            tube_idx * num_detectors_in_tube:
                num_detectors_in_tube + tube_idx*num_detectors_in_tube]

        # check how many detectors are masked in the tube
        # same size as the tube, True where masked
        detectors_masked = np.array(
            [detector_info_input_ws.isMasked(detector_id_to_index[id])
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

                output_ws.setY(ws_index_to_edit, np.array([y_new_value]))
                output_ws.setE(ws_index_to_edit, np.array([e_new_value]))
                detector_info_output_ws.setMasked(ws_index_to_edit, False)

            # mask pixels in WS Out where detectors_inf
            for detector_id in tube_detector_ids[detectors_inf]:
                ws_index_to_edit = detector_id_to_index[detector_id]
                detector_info_output_ws.setMasked(ws_index_to_edit, True)
                output_ws.setY(ws_index_to_edit, np.array([1]))
                output_ws.setE(ws_index_to_edit, np.array([0]))

        elif num_of_detectors_masked > min_detectors_per_tube:
            logger.error("Skipping tube {}. Too many masked or dead pixels."
                         "".format(tube_idx))
        # Another iteration
        ws_index += num_detectors_in_tube
    return output_ws


def _get_detector_id_to_ws_index(flood_ws):
    spectrum_info = flood_ws.spectrumInfo()
    detector_id_to_index = OrderedDict(
        (flood_ws.getSpectrum(ws_index).getDetectorIDs()[0], ws_index)
        for ws_index in range(flood_ws.getNumberHistograms())
        if not spectrum_info.isMonitor(ws_index))
    return detector_id_to_index


def inf_value_to_mask(ws):
    """In case the input workspace has infinite values (EMPTY_DBL)
    Mask these pixels and set Value to 1 and Error to 0.
    """

    assert ws.getNumberBins() == 1, "This only supports integrated WS"

    detector_info_input_ws = ws.detectorInfo()
    for i in range(ws.getNumberHistograms()):
        if ws.readY(i)[0] == Property.EMPTY_DBL:
            detector_info_input_ws.setMasked(i, True)
            ws.setY(i, np.array([1]))
            ws.setE(i, np.array([0]))
