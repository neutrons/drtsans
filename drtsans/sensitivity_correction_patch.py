import numpy as np
import os
from mantid.kernel import Property, logger
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.settings import unique_workspace_name as uwn
r"""
Links to mantid algorithms
https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
https://docs.mantidproject.org/nightly/algorithms/SaveNexusProcessed-v1.html
https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html
https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
"""
from mantid.simpleapi import mtd, CloneWorkspace, DeleteWorkspace, SaveNexusProcessed, Integration, CreateWorkspace

from drtsans import detector

__all__ = ['calculate_sensitivity_correction']


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


def calculate_sensitivity_correction(input_workspace, min_threshold=0.5, max_threshold=2.0,
                                     poly_order=2, min_detectors_per_tube=50, filename=None,
                                     component_name='detector1', output_workspace=None):
    """
    Calculate the detector sensitivity

    Prerequisites for input workspace:
    1. All previously masked values to NaN as required by Numpy functions but not masked pixels for beam centers
    2. All masked pixels at beam centers are set to -infinity

    **Mantid algorithms used:**
    :ref:`SaveNexusProcessed <algm-SaveNexusProcessed-v1>`
    :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace to calculate the sensitivity from
    min_threshold: float
        Minimum threshold for efficiency value.
    max_threshold: float
        Maximum threshold for efficiency value
    poly_order : int
        ploy order.  default to 2
    min_detectors_per_tube : int, optional
        Minimum detectors with a value existing in the tube to fit. Only fits
        tubes with at least `min_detectors_per_tube` (the default is 50).
    component_name : str, optional
        Component name to  (the default is 'detector1', which is the main
        detector)
    filename: str
        Name of the file to save the sensitivity calculation to
    output_workspace: ~mantid.api.MatrixWorkspace
        The calculated sensitivity workspace
    """
    if output_workspace is None:
        output_workspace = '{}_sensitivity'.format(input_workspace)

    # Wavelength bins are summed together to remove the time-of-flight nature according
    # to equations A3.1 and A3.2
    input_workspace = mtd[str(input_workspace)]

    if input_workspace.blocksize() != 1:
        input_workspace = Integration(InputWorkspace=input_workspace, OutputWorkspace=uwd())
        delete_input_workspace = True
    else:
        delete_input_workspace = False

    # The average and uncertainty in the average are determined from the masked pattern
    # according to equations A3.3 and A3.4
    # numpy.flatten() used to more easily find the mean and uncertainty using numpy.
    y = input_workspace.extractY().flatten()
    y_uncertainty = input_workspace.extractE().flatten()

    n_elements =\
        input_workspace.getNumberHistograms() - np.count_nonzero(np.isnan(y)) - np.count_nonzero(np.isneginf(y))
    F = np.sum([value for value in y if not np.isnan(value) and not np.isneginf(value)])/n_elements
    dF = np.sqrt(np.sum([value**2 for value in y_uncertainty
                         if not np.isnan(value) and not np.isneginf(value)]))/n_elements
    II = y/F
    dI = II * np.sqrt(np.square(y_uncertainty/y) + np.square(dF/F))

    # Any pixel in II less than min_threshold or greater than max_threshold is masked
    for i, value in enumerate(II):
        if not np.isnan(value) and not np.isneginf(value):
            if value < min_threshold or value > max_threshold:
                II[i] = np.nan
                dI[i] = np.nan

    comp = detector.Component(input_workspace, component_name)  # 'detector` except wing detector

    # The next step is to fit the data in each tube with a second order polynomial as shown in
    # Equations A3.9 and A3.10. Use result to fill in NaN values.
    num_interpolated_tubes = 0
    for j in range(0, comp.dim_y):
        xx = []
        yy = []
        ee = []
        # beam center masked pixels
        masked_indices = []
        for i in range(0, comp.dim_x):
            index = comp.dim_x*j + i
            if np.isneginf(II[index]):
                masked_indices.append([i, index])
            elif not np.isnan(II[index]):
                xx.append(i)
                yy.append(II[index])
                ee.append(dI[index])

        if len(masked_indices) == 0:
            # no masked/centermasked pixels
            # no need to do interpolation
            continue
        # This shall be an option later
        if len(xx) < min_detectors_per_tube:
            logger.error("Skipping tube with indices {} with {} non-masked value. Too many "
                         "masked or dead pixels.".format(j, len(xx)))
            print('Tube {} .................... Valid Pixels = {}......... Skip'.format(j, len(xx)))
            continue

        # Do poly fit
        num_interpolated_tubes += 1
        polynomial_coeffs, cov_matrix = np.polyfit(xx, yy, poly_order, w=np.array(ee), cov=True)

        # Errors in the least squares is the sqrt of the covariance matrix
        # (correlation between the coefficients)
        e_coeffs = np.sqrt(np.diag(cov_matrix))

        masked_indices = np.array(masked_indices)

        # The patch is applied to the results of the previous step to produce S2(m,n).
        y_new = np.polyval(polynomial_coeffs, masked_indices[:, 0])
        # errors of the polynomial
        e_new = np.sqrt(e_coeffs[2]**2 + (e_coeffs[1]*masked_indices[:, 0])**2 +
                        (e_coeffs[0]*masked_indices[:, 0]**2)**2)
        for i, index in enumerate(masked_indices[:, 1]):
            II[index] = y_new[i]
            dI[index] = e_new[i]

    # The final sensitivity, S(m,n), is produced by dividing this result by the average value
    # per Equations A3.13 and A3.14
    # numpy.flatten() used to more easily find the mean and uncertainty using numpy.
    n_elements =\
        input_workspace.getNumberHistograms() - np.count_nonzero(np.isnan(II)) - np.count_nonzero(np.isneginf(II))
    F = np.sum([value for value in II if not np.isnan(value) and not np.isneginf(value)])/n_elements
    dF = np.sqrt(np.sum([value**2 for value in dI if not np.isnan(value) and not np.isneginf(value)]))/n_elements
    output = II/F
    output_uncertainty = output * np.sqrt(np.square(dI/II) + np.square(dF/F))

    CreateWorkspace(DataX=[1., 2.],
                    DataY=output,
                    DataE=output_uncertainty,
                    Nspec=input_workspace.getNumberHistograms(),
                    UnitX='wavelength',
                    OutputWorkspace=output_workspace)
    if delete_input_workspace:
        DeleteWorkspace(input_workspace)
    if filename:
        path = os.path.join(os.path.expanduser("~"), filename)
        SaveNexusProcessed(InputWorkspace=output_workspace, Filename=path)
    return mtd[output_workspace]
