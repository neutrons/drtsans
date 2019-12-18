from collections import OrderedDict

import numpy as np
import os

from mantid.kernel import Property, logger
from mantid.simpleapi import mtd, CloneWorkspace, CalculateEfficiency,\
    DeleteWorkspace, Divide, LoadNexusProcessed, MaskDetectors, \
    MaskDetectorsIf, SaveNexusProcessed, CreateSingleValuedWorkspace, Quadratic, \
    Fit
from drtsans.path import exists as path_exists
from drtsans.settings import unique_workspace_name as uwd
from drtsans import detector

__all__ = ['apply_sensitivity_correction', 'calculate_sensitivity_correction', 'prepare_sensitivity_correction']


class Detector(object):
    r"""
    Auxiliary class that has all the information about a detector
    It allows to read tube by tube.
    """

    def __init__(self, workspace, component_name):
        self._workspace = workspace
        self._current_start_ws_index = None  # first workspace index of the currently considered tube
        self._current_stop_ws_index = None  # last workspace index of the currently considered tube
        self._tube_ws_indices = None  # iterator for tube endpoints, given as workspace indexes

        # Public variables
        self.component_name = None  # name of the assembly of pixels making up the Detector
        self.n_tubes = None  # number of tubes in the detector
        self.n_pixels_per_tube = None
        self.first_det_id = None  # pixel ID for first pixel detector that is not a monitor
        self.last_det_id = None  # pixel ID for last pixel detector
        self.detector_id_to_ws_index = None # mapping from pixel ID to workspace index
        self.data_y = None
        self.data_e = None
        self.tube_ws_indices = None

        # Initialize attributes of the Detector
        self._detector_details(component_name)
        self._detector_id_to_ws_index()
        self._extract_data()
        self._set_tube_ws_indices()

    def _detector_details(self, component_name):
        r"""
        Initializes the following attributes of the Detector:
        component_name, n_tubes, n_pixels_per_tube, first_det_id, last_det_id

        Parameters
        ----------
        component_name : string
            Name of the pixel assembly, one component of the instrument.
        """
        self.component_name = component_name
        i = self._workspace.getInstrument()
        component = i.getComponentByName(component_name)
        num_pixels = 1
        # dive into subelements until get a detector
        while component.type() != 'DetectorComponent' and component.type() != 'GridDetectorPixel':
            self.n_pixels_per_tube = component.nelements()
            num_pixels *= self.n_pixels_per_tube
            component = component[0]
        self.first_det_id = component.getID()
        self.last_det_id = self.first_det_id + num_pixels - 1
        self.n_tubes = int(num_pixels/self.n_pixels_per_tube)

    def _detector_id_to_ws_index(self):
        r"""
        Maps the detector ID of one pixel to one workspace index.
        Initializes attribute ``detector_id_to_ws_index``.
        """

        spectrum_info = self._workspace.spectrumInfo()
        detector_id_to_index = []
        for ws_index in range(self._workspace.getNumberHistograms()):
            if spectrum_info.isMonitor(ws_index) is True:
                continue
            detector_id_to_index.append((self._workspace.getSpectrum(ws_index).getDetectorIDs()[0], ws_index))
        self.detector_id_to_ws_index = OrderedDict(detector_id_to_index)

    def _extract_data(self):
        r"""
        Extract intensitites and associated uncertainties from the workspace.
        Initializes attributes ``data_y`` and ``data_e``
        """
        self.data_y = self._workspace.extractY()
        self.data_e = self._workspace.extractE()

    def _set_tube_ws_indices(self):
        r"""
        Initializes attribute ``_tube_ws_indices`` by assigning to it an iterator that yields pairs of the form
        ``(start, end)``, where ``start`` and ``end`` are the starting and ending workspace indexes for a given tube.
        The iterator yields as many pairs as tubes in Detector.
        """
        tube_ws_indices = []
        for tube_idx in range(self.n_tubes):
            first_det_id = self.first_det_id + tube_idx * self.n_pixels_per_tube  # tube starts with this pixel ID
            last_det_id = first_det_id + self.n_pixels_per_tube - 1  # tube ends with this pixel ID

            first_ws_index = self.detector_id_to_ws_index[first_det_id]  # tube starts with this workspace index
            last_ws_index = self.detector_id_to_ws_index[last_det_id]  # tube ends with this workspace index

            tube_ws_indices.append((first_ws_index, last_ws_index))
        self._tube_ws_indices = iter(tube_ws_indices)  # return an iterator

    def next_tube(self):
        r"""Initializes/ updates attributes ``_current_start_ws_index`` and ``_current_stop_ws_index``"""
        self._current_start_ws_index, self._current_stop_ws_index = next(self._tube_ws_indices)

    def get_current_ws_indices(self):
        r"""
        First and last workspace indices for the currently considered tube.

        Returns
        -------
        tuple
        """
        return self._current_start_ws_index, self._current_stop_ws_index

    def get_current_ws_indices_range(self):
        r"""
        Array of workspace indexes for the currently considered tube.

        Returns
        -------
        ~numpy.ndarray
        """

        return np.array(range(self._current_start_ws_index, self._current_stop_ws_index+1))

    def get_ws_data(self):
        r"""
        Intensities and associated uncertainties for the currently considered tube.

        Returns
        -------
        tuple
            A two-item tuple containing, in this order, intensites and uncertainties in the shape of ~numpy.ndarray.
        """
        return (self.data_y[self._current_start_ws_index: self._current_stop_ws_index+1].flatten(),
                self.data_e[self._current_start_ws_index: self._current_stop_ws_index+1].flatten())

    def get_pixels_masked(self):
        r"""
        Pixel masks for the currently considered tube.

        Returns
        -------
        ~numpy.ndarray
            Array of ``Bool`` values, with :py:obj:`True` for the masked pixels and :py:obj:`False` otherwise.
        """
        spectrum_info = self._workspace.spectrumInfo()
        return np.array([spectrum_info.isMasked(int(idx)) for idx in self.get_current_ws_indices_range()])

    def get_pixels_infinite(self):
        r"""
        Pixel mask for pixels with non-finite intensities in the currently considered tube.
        Returns an array of booleans for this tube
        where the pixel count is EMPTY_DBL

        Returns
        -------
        ~numpy.ndarray
            Array of ``Bool`` values, with :py:obj:`True` for the pixels with non-finite intensities, and
            :py:obj:`False` otherwise.
        """
        return np.array([self._workspace.readY(int(idx))[0] == Property.EMPTY_DBL
                         for idx in self.get_current_ws_indices_range()])

    def get_y_coordinates(self):
        r"""
        Y-coordinates of the pixels for the currently considered tube.

        Returns
        -------
        ~numpy.ndarray
        """
        detector_info = self._workspace.spectrumInfo()
        return np.array([detector_info.position(int(idx))[1] for idx in self.get_current_ws_indices_range()])


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
        flood_ws, OutputWorkspace=uwn(
            prefix="__sensitivity_"))
    detector_info_output_ws = output_ws.spectrumInfo()

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

    detector_info_input_ws = ws.spectrumInfo()
    for i in range(ws.getNumberHistograms()):
        if ws.readY(i)[0] == Property.EMPTY_DBL:
            detector_info_input_ws.setMasked(i, True)
            ws.setY(i, np.array([1]))
            ws.setE(i, np.array([0]))


# flake8: noqa: C901
def apply_sensitivity_correction(input_workspace, sensitivity_filename=None,
                                 sensitivity_workspace=None,
                                 min_threshold=None,
                                 max_threshold=None,
                                 output_workspace=None):
    '''Apply a previously calculated sensitivity correction

    **Mantid algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,
    :ref:`Divide <algm-Divide-v1>`,
    :ref:`LoadNexusProcessed <algm-LoadNexusProcessed-v1>`,
    :ref:`MaskDetectors <algm-MaskDetectors-v1>`
    :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        workspace to apply the correction to
    sensitivity_filename: str
        file containing previously calculated sensitivity correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    min_threshold: float or None
        if not None, the data will be masked if the sensitivity
        is below this threshold
    max_threshold: float or None
        if not None, the data will be masked if the sensitivity
        is above this threshold
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

    # additional masking dependent on threshold
    temp_sensitivity = CloneWorkspace(InputWorkspace=sensitivity_workspace,
                                      OutputWorkspace=uwn(prefix="__sensitivity_"))
    if min_threshold is not None:
        MaskDetectorsIf(InputWorkspace=temp_sensitivity,
                        Operator='LessEqual',
                        Value=min_threshold,
                        OutputWorkspace=temp_sensitivity)
    if max_threshold is not None:
        MaskDetectorsIf(InputWorkspace=temp_sensitivity,
                        Operator='GreaterEqual',
                        Value=max_threshold,
                        OutputWorkspace=temp_sensitivity)
    Divide(LHSWorkspace=output_workspace, RHSWorkspace=temp_sensitivity,
           OutputWorkspace=output_workspace)
    DeleteWorkspace(temp_sensitivity)

    if cleanupSensitivity:
        DeleteWorkspace(sensitivity_workspace)

    return mtd[output_workspace]


def calculate_sensitivity_correction(input_workspace, min_threashold=0.5, max_threshold=2.0,
                                     filename=None, output_workspace=None):
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
        Minimum threshold for efficiency value.
    max_threashold: float
        Maximum threshold for efficiency value
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    '''
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    CalculateEfficiency(InputWorkspace=input_workspace, OutputWorkspace=output_workspace,
                        MinThreshold=min_threashold, MaxThreshold=max_threshold)
    MaskDetectorsIf(InputWorkspace=output_workspace, OutputWorkspace=output_workspace,
                    Mode='SelectIf', Operator='Equal', Value=Property.EMPTY_DBL)
    if filename is not None:
        SaveNexusProcessed(InputWorkspace=output_workspace, Filename=filename)
    return mtd[output_workspace]

def prepare_sensitivity_correction(input_workspace,  min_threshold=0.5,  max_threshold=2.0,
                                     filename=None,  output_workspace=None):
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
        Minimum threshold for efficiency value.
    max_threashold: float
        Maximum threshold for efficiency value
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    '''
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    y = input_workspace.extractY().flatten()
    indices_to_mask = np.arange(len(y))[np.isnan(y)]
    original_mask = np.isnan(y)
    F = np.nanmean(y)
    MaskDetectors(input_workspace, WorkspaceIndexList=indices_to_mask)
    n_elements = 0
    for i in range(input_workspace.getNumberHistograms()):
        n_elements += len(input_workspace.readY(i))
    n_elements -= len(indices_to_mask)
    y_uncertainty = input_workspace.extractE().flatten()
    dF = np.sqrt(np.nansum(np.power(y_uncertainty, 2)))/n_elements
    F_ws = CreateSingleValuedWorkspace(DataValue=F, ErrorValue=dF, OutputWorkspace=uwd())
    II = Divide(LHSWorkspace=input_workspace, RHSWorkspace=F_ws, OutputWorkspace=uwd())

    MaskDetectorsIf(InputWorkspace=II, OutputWorkspace=II,
                    Mode='SelectIf', Operator='Greater', Value=max_threshold)

    MaskDetectorsIf(InputWorkspace=II, OutputWorkspace=II,
                    Mode='SelectIf', Operator='Less', Value=min_threshold)

    det_info = II.detectorInfo()
    comp = detector.Component(II, 'detector1')

    for j in range(0, comp.dim_y):
        xx = []
        yy = []
        ee = []
        masked_indices = []
        for i in range(0, comp.dim_x):
            index = comp.dim_x*j + i
            if det_info.isMasked(index):
                masked_indices.append([i, index])
            else:
                xx.append(i)
                yy.append(II.readY(index)[0])
                ee.append(II.readE(index)[0])

        polynomial_coeffs, cov_matrix = np.polyfit(xx, yy, 2, w=np.array(ee), cov=True)
        # Errors in the least squares is the sqrt of the covariance matrix
        # (correlation between the coefficients)
        e_coeffs = np.sqrt(np.diag(cov_matrix))
        masked_indices = np.array(masked_indices)

        y_new = np.polyval(polynomial_coeffs, masked_indices[:, 0])
        a = masked_indices[:, 0]
        # errors of the polynomial
        e_new = np.sqrt(e_coeffs[2]**2 + (e_coeffs[1]*masked_indices[:, 0])**2 +
                        (e_coeffs[0]*masked_indices[:, 0]**2)**2)
        for i, index in enumerate(masked_indices[:,1]):
            if original_mask[index]:
              det_info.setMasked(int(index), False)
              II.setY(int(index), [y_new[i]])
              II.setE(int(index), np.array(e_new[i]))
            else:
              II.setY(int(index), [np.nan])
              II.setE(int(index), np.array(np.nan))

    y = II.extractY().flatten()
    indices_to_mask = np.arange(len(y))[np.isnan(y)]
    F = np.nanmean(y)
    n_elements = 0
    for i in range(input_workspace.getNumberHistograms()):
        n_elements += len(input_workspace.readY(i))
    n_elements -= len(indices_to_mask)
    y_uncertainty = II.extractE().flatten()
    dF = np.sqrt(np.nansum(np.power(y_uncertainty, 2)))/n_elements
    F_ws = CreateSingleValuedWorkspace(DataValue=F, ErrorValue=dF, OutputWorkspace=uwd())
    output_workspace = Divide(LHSWorkspace=II, RHSWorkspace=F_ws, OutputWorkspace=uwd())
    return output_workspace