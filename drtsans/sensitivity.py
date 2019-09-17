from __future__ import absolute_import, division, print_function

from collections import OrderedDict

import numpy as np
import os

from mantid.kernel import Property, logger
from mantid.simpleapi import mtd, CloneWorkspace, CalculateEfficiency,\
    DeleteWorkspace, Divide, LoadNexusProcessed, MaskDetectors, \
    MaskDetectorsIf, SaveNexusProcessed
from ornl.path import exists as path_exists
from ornl.settings import unique_workspace_name


class Detector(object):
    '''
    Auxiliary class that has all the information about a detector
    It allows to read tube by tube.
    '''

    def __init__(self, workspace, component_name):
        self._workspace = workspace
        self._current_start_ws_index = None
        self._current_stop_ws_index = None
        self._tube_y_coordinates = None

        # Public variables
        self.n_tubes = None
        self.n_pixels_per_tube = None
        self.first_det_id = None
        self.last_det_id = None
        self.data_y = None
        self.data_e = None
        self.tube_ws_indices = None

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
        '''Sets `tube_ws_indices` as an iterator for the next round of tube
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
        '''For the current tube returns the first and last workspace indices

        Returns
        -------
        tuple
        '''
        return (self._current_start_ws_index, self._current_stop_ws_index)

    def get_current_ws_indices_range(self):
        '''For the current tube returns the range delimited by the first and
        last workspace indices.

        Returns
        -------
        numpy.array
        '''

        return np.array(range(self._current_start_ws_index,
                              self._current_stop_ws_index))

    def get_ws_data(self):
        '''Returns the current tube data and error

        Returns
        -------
        Tuple of arrays
            Y and Error
        '''
        return (self.data_y[self._current_start_ws_index:
                            self._current_stop_ws_index].flatten(),
                self.data_e[self._current_start_ws_index:
                            self._current_stop_ws_index].flatten())

    def get_pixels_masked(self):
        '''Returns an array of booleans for this tube
        where the detector is masked

        Returns
        -------
        np.array
        '''
        detector_info = self._workspace.detectorInfo()
        return np.array([detector_info.isMasked(idx) for idx in range(
            self._current_start_ws_index, self._current_stop_ws_index)])

    def get_pixels_infinite(self):
        '''Returns an array of booleans for this tube
        where the pixel count is EMPTY_DBL

        Returns
        -------
        np.array
        '''
        return np.array([self._workspace.readY(idx)[0] == Property.EMPTY_DBL
                         for idx in range(self._current_start_ws_index,
                                          self._current_stop_ws_index)])

    def get_y_coordinates(self):
        '''Return a numpy array of the current tube Y coordinates
        The Y coordinates are allways the same no matter the position of the
        tube, so cache it `self._tube_y_coordinates`

        Returns
        -------
        np.array
        '''
        if self._tube_y_coordinates is None:
            detector_info = self._workspace.detectorInfo()
            self._tube_y_coordinates = np.array([
                detector_info.position(idx)[1] for idx in range(
                    self._current_start_ws_index, self._current_stop_ws_index)
            ])
        return self._tube_y_coordinates


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
    assert flood_ws.blocksize() == 1, "This only supports integrated WS"

    d = Detector(flood_ws, component_name)
    # Lets get the output workspace
    output_ws = CloneWorkspace(
        flood_ws, OutputWorkspace=unique_workspace_name(
            prefix="__sensitivity_"))
    detector_info_output_ws = output_ws.detectorInfo()

    while True:
        try:
            d.next_tube()
        except StopIteration:
            break  # Iterator exhausted: stop the loop

        x = d.get_y_coordinates()
        detectors_masked = d.get_pixels_masked()
        detectors_inf = d.get_pixels_infinite()
        # Count the detectors masked or infinite (True)
        num_of_detectors_masked = np.count_nonzero(detectors_masked)
        num_of_detectors_inf = np.count_nonzero(detectors_inf)

        if num_of_detectors_masked > 0 and \
            (d.n_pixels_per_tube - num_of_detectors_masked -
                num_of_detectors_inf) > min_detectors_per_tube:
            # number of detectors with values > min_detectors_per_tube

            # Let's fit
            y, e = d.get_ws_data()
            y_new, e_new = _interpolate_tube(
                x, y, e, detectors_masked, detectors_inf, polynomial_degree)

            # Set output workspace with interpolated values
            for ws_index, y_new_value, e_new_value in zip(
                    d.get_current_ws_indices_range()[detectors_masked],
                    y_new[detectors_masked], e_new[detectors_masked]):

                output_ws.setY(int(ws_index), np.array([y_new_value]))
                output_ws.setE(int(ws_index), np.array([e_new_value]))
                detector_info_output_ws.setMasked(int(ws_index), False)

            # mask pixels in WS Out where detectors_inf
            for ws_index in d.get_current_ws_indices_range()[detectors_inf]:
                detector_info_output_ws.setMasked(int(ws_index), True)
                output_ws.setY(int(ws_index), np.array([1]))
                output_ws.setE(int(ws_index), np.array([0]))

        elif num_of_detectors_masked > min_detectors_per_tube:
            logger.error("Skipping tube with indices {}. Too many "
                         "masked or dead pixels.".format(
                             d.get_current_ws_indices()))
    return output_ws


def inf_value_to_mask(ws):
    """In case the input workspace has infinite values (EMPTY_DBL)
    Mask these pixels and set Value to 1 and Error to 0.
    """

    assert ws.blocksize() == 1, "This only supports integrated WS"

    detector_info_input_ws = ws.detectorInfo()
    for i in range(ws.getNumberHistograms()):
        if ws.readY(i)[0] == Property.EMPTY_DBL:
            detector_info_input_ws.setMasked(i, True)
            ws.setY(i, np.array([1]))
            ws.setE(i, np.array([0]))


def apply_sensitivity_correction(input_workspace, sensitivity_filename=None,
                                 sensitivity_workspace=None,
                                 output_workspace=None):
    '''Apply a previously calculated sensitivity correction

    **Mantid algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,
    :ref:`Divide <algm-Divide-v1>`,
    :ref:`LoadNexusProcessed <algm-LoadNexusProcessed-v1>`,
    :ref:`MaskDetectors <algm-MaskDetectors-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        workspace to apply the correction to
    sensitivity_filename: str
        file containing previously calculated sensitivity correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    output_workspace:  ~mantid.api.MatrixWorkspace
        corrected workspace. This is the input workspace by default
    '''
    if output_workspace is None:
        output_workspace = str(input_workspace)

    cleanupSensitivity = False
    if sensitivity_workspace is None or \
            str(sensitivity_workspace) not in mtd:  # load the file
        if sensitivity_workspace is None:
            sensitivity_workspace = os.path.split(sensitivity_filename)[-1]
            sensitivity_workspace = sensitivity_workspace.split('.')[0]
        cleanupSensitivity = True
        if not path_exists(sensitivity_filename):
            msg = 'Cannot find file "{}"'.format(sensitivity_filename)
            raise RuntimeError(msg)
        LoadNexusProcessed(Filename=sensitivity_filename,
                           OutputWorkspace=sensitivity_workspace)
    if (not sensitivity_workspace) or (str(sensitivity_workspace) not in mtd):
        raise RuntimeError('No sensitivity workspace provided')

    if str(input_workspace) != str(output_workspace):
        CloneWorkspace(InputWorkspace=input_workspace,
                       OutputWorkspace=output_workspace)
    MaskDetectors(Workspace=output_workspace,
                  MaskedWorkspace=sensitivity_workspace)
    Divide(LHSWorkspace=output_workspace, RHSWorkspace=sensitivity_workspace,
           OutputWorkspace=output_workspace)

    if cleanupSensitivity:
        DeleteWorkspace(sensitivity_workspace)

    return mtd[output_workspace]


def calculate_sensitivity_correction(input_workspace, min_threashold,
                                     max_threshold, filename,
                                     output_workspace=None):
    '''
    Calculate the detector sensitivity

    **Mantid algorithms used:**
    :ref:`CalculateEfficiency <algm-CalculateEfficiency-v1>`,
    :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
    :ref:`SaveNexusProcessed <algm-SaveNexusProcessed-v1>`


    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace to calculate the sensitivity from
    min_threashold: float
        Minimum threshold for sensitivity value
    max_threashold: float
        Maximum threshold for sensitivity value
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    '''
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    CalculateEfficiency(InputWorkspace=str(input_workspace),
                        OutputWorkspace=output_workspace,
                        MinThreshold=min_threashold,
                        MaxThreshold=max_threshold)
    MaskDetectorsIf(InputWorkspace=output_workspace,
                    OutputWorkspace=output_workspace, Mode='SelectIf',
                    Operator='Equal', Value=Property.EMPTY_DBL)
    SaveNexusProcessed(InputWorkspace=output_workspace,
                       Filename=filename)
    return mtd[output_workspace]
