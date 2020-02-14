# This module contains legacy code to do sensitivity correction with detector patching method.
# It is not supposed to be used or supported after the new algorithm passes acceptance.
import numpy as np
r"""
Links to mantid algorithms
https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html
https://docs.mantidproject.org/nightly/algorithms/SaveNexusProcessed-v1.html
"""
from mantid.simpleapi import mtd, CalculateEfficiency, MaskDetectorsIf, SaveNexusProcessed, CloneWorkspace
from mantid.kernel import Property
from drtsans.settings import unique_workspace_name as uwn
from drtsans.detector import Detector
# from mantid.simpleapi import mtd,  DeleteWorkspace, SaveNexusProcessed, Integration, CreateWorkspace

__all__ = ['calculate_sensitivity_correction']


# This is old sensitivity calculation method.  It is kept temporarily for comparison purpose
def calculate_sensitivity_correction(input_workspace, min_threshold=0.5, max_threshold=2.0,
                                     filename=None, output_workspace=None):
    """
    Calculate the detector sensitivity

    **Mantid algorithms used:**
    :ref:`CalculateEfficiency <algm-CalculateEfficiency-v1>`,
    :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
    :ref:`SaveNexusProcessed <algm-SaveNexusProcessed-v1>`


    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace to calculate the sensitivity from
    min_threshold: float
        Minimum threshold for efficiency value.
    max_threshold: float
        Maximum threshold for efficiency value
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    """
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    CalculateEfficiency(InputWorkspace=input_workspace, OutputWorkspace=output_workspace,
                        MinThreshold=min_threshold, MaxThreshold=max_threshold)
    MaskDetectorsIf(InputWorkspace=output_workspace, OutputWorkspace=output_workspace,
                    Mode='SelectIf', Operator='Equal', Value=Property.EMPTY_DBL)
    if filename is not None:
        SaveNexusProcessed(InputWorkspace=output_workspace, Filename=filename)
    return mtd[output_workspace]


# The following methods are only called by several test_sensitivity.py
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
