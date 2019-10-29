import numpy as np
from mantid.api import mtd

r"""
Hyperlinks to Mantid algorithms
CloneWorkspace <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
Fit <https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html>
Plus <https://docs.mantidproject.org/nightly/algorithms/Plus-v1.html>
SaveNexus <https://docs.mantidproject.org/nightly/algorithms/SaveNexus-v1.html>
"""
from mantid.simpleapi import CloneWorkspace, Fit, Plus, SaveNexus
r"""
Hyperlinks to drtsans functions
namedtuplefy, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
calculate_transmission <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/transmission.py>
clipped_bands_from_logs, transmitted_bands available at: 
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans
/correct_frame.py>
sample_aperture_diameter, source_aperture_diameter available at:
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/geometry.py>
"""  # noqa: E501
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.transmission import calculate_transmission as calculate_raw_transmission
from drtsans.tof.eqsans.correct_frame import clipped_bands_from_logs, transmitted_bands
from drtsans.tof.eqsans.geometry import sample_aperture_diameter, source_aperture_diameter
from drtsans.geometry import sample_detector_distance, source_sample_distance


# Symbols to be exported to the eqsans namespace
__all__ = ['beam_radius', 'calculate_transmission', 'fit_raw_transmission']


def calculate_transmission(input_sample, input_reference, radius=None, radius_unit='mm',
                           fit_function='name=UserFunction,Formula=a*x+b', output_workspace=None):
    """
    Calculate the transmission coefficients at zero scattering angle
    from already prepared sample and reference data.

    For EQ-SANS, one additional step fitting the raw values with an
    analytic wavelength-dependent function may be requested

    Parameters
    ----------
    input_sample: str, ~mantid.api.MatrixWorkspace
        Prepared sample workspace (possibly obtained with an attenuated beam)
    input_reference: str, ~mantid.api.MatrixWorkspace
        Prepared direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        Radius around the bean center for pixel integration, in milimeters.
        If None, radius will be obtained or calculated using ``input_reference``.
    radius_unit: str
        Either 'mm' or 'm', and only used in conjunction with option ``radius``.
    fit_function: str
        String representation of the fit function. See Mantid's
        :ref:`UserFunction <func-UserFunction>` or any of Mantid's built-in functions.
        The default value represents a linear model. If this option is left as :py:obj:`None` or empty
        string, then no fitting is performed.
    output_workspace: str
        Name of the output workspace containing the raw transmission values.
        If :py:obj:`None`, an anonymous hidden name will be automatically provided.

    Returns
    -------
    MatrixWorkspace
        Workspace containing the transmission values
    """
    if output_workspace is None:
        output_workspace = unique_workspace_dundername()
    transmission_workspace = calculate_raw_transmission(input_sample, input_reference, radius=radius,
                                                        radius_unit=radius_unit, output_workspace=output_workspace)
    if bool(fit_function) is True:
        transmission_workspace = fit_raw_transmission(transmission_workspace, fit_function=fit_function).transmission
    return transmission_workspace


def beam_radius(input_workspace, unit='mm'):
    r"""
    Calculate the beam radius impinging on the detector bank.

    .. math::

           R_{beam} = R_{sampleAp} + SDD * (R_{sampleAp} + R_{sourceAp}) / SSD

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace, str
        Input workspace, contains all necessary info in the logs
    unit: str
        Either 'mm' or 'm'

    Returns
    -------
    float
        Estimated beam radius
    """
    sample_aperture_diam = sample_aperture_diameter(input_workspace, unit=unit)
    source_aperture_diam = source_aperture_diameter(input_workspace, unit=unit)
    source_sample_dist = source_sample_distance(input_workspace)
    sample_detector_dist = sample_detector_distance(input_workspace)
    return sample_aperture_diam +\
        sample_detector_dist * (sample_aperture_diam + source_aperture_diam) / (2 * source_sample_dist)


@namedtuplefy
def fit_band(input_workspace, band, fit_function='name=UserFunction,Formula=a*x+b',
             output_workspace=None):
    r"""
    Fit the wavelength dependence of the raw zero-angle transmission
    values with a function within a wavelength band.

    Parameters
    ----------
    input_workspace: MatrixWorkspace, str
        Input workspace containing the raw transmission values
    band: Wband
        Wavelength band over which to carry out the fit
    fit_function: str
        String representation of the fit function. See Mantid's
        `UserFunction` or any of Mantid's fit functions
    output_workspace: str
        Name of the output workspace containing fitted transmission
        values within the band, and zero elsewhere. If `None`, a random
        hidden name will be generated for the output workspace.

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - fitted_workspace: MatrixWorkspace, transmission values within the band,
          zero elsewhere. The name of this workspace is `output_workspace`
        - mantid_fit_output: namedtuple, output of calling Mantid's Fit algorithm
    """
    if output_workspace is None:
        output_workspace = unique_workspace_dundername()

    ws = mtd[str(input_workspace)]

    # Carry out the fit only over the wavelength band
    x, y = ws.dataX(0), ws.dataY(0)
    band_indexes = np.where((x >= band.min) & (x < band.max))[0]

    # It has been found that raw transmission values at the edges of the wavelength band do fall short of the
    # expected values. These anomalous values should be excluded from the fit
    min_y = 1e-3 * np.mean(y[band_indexes[:-1]])  # 1e-3 pure heuristics
    # Find the low-wavelength edge with non-zero intensities
    i = 0
    while x[i] < band.min or y[i] < min_y:
        i += 1
    lower_bin_boundary = i
    # Find the high-wavelength edge with non-zero intensities
    while i < len(y) and x[i] < band.max and y[i] > min_y:
        i += 1
    upper_bin_boundary = i
    start_x, end_x = x[lower_bin_boundary], x[upper_bin_boundary - 1]  # this is our restricted fitting range

    mantid_fit = Fit(Function=fit_function, InputWorkspace=ws.name(), WorkspaceIndex=0,
                     StartX=start_x, EndX=end_x, Output=unique_workspace_dundername())

    # Insert the fitted band into the wavelength range of the input workspace
    fitted = CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)
    insert_fitted_values(output_workspace, mantid_fit, lower_bin_boundary, upper_bin_boundary)
    insert_fitted_errors(output_workspace, mantid_fit, lower_bin_boundary, upper_bin_boundary)
    return dict(fitted_workspace=fitted, mantid_fit_output=mantid_fit)


@namedtuplefy
def fit_raw_transmission(input_workspace, fit_function='name=UserFunction,Formula=a*x+b',
                         output_workspace=None):
    r"""
    Fit the wavelength dependence of the raw zero-angle transmission
    values with a model.

    If working in frame skipping mode, apply the fit separately to the
    wavelength bands of the lead and skipped pulses.

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Input workspace containing the raw transmission values
    output_workspace: str
        Name of the workspace containing the fitted transmission values
        and errors. If None, the input worskpace will be overwritten
    fit_function: str
        String representation of the fit function. See Mantid's
        `UserFunction` or any of Mantid's fit functions

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - transmission: workspace containing the fitted transmissio values and errors. The name of this workspace
        is the value of `output_workspace`
        - lead_transmission: workspace containing the fitted transmission values and
            errors of the lead pulse
        - lead_mantid_fit: return value of running Mantid's Fit algorithm when
            fitting the raw transmission over the lead pulse wavelength range
        - skip_transmission: workspace containing the fitted transmission values and
            errors of the skip pulse. None if not working in frame skipping
            mode
        - skip_mantid_fit: return value of running Mantid's Fit algorithm when
            fitting the raw transmission over the skip pulse wavelength range
            None f not working in frame skipping mode
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    # Fit only over the range of the transmitted wavelength band(s)
    try:
        wavelength_bands = clipped_bands_from_logs(input_workspace)
    except AttributeError:  # bands are not logged
        wavelength_bands = transmitted_bands(input_workspace)  # calculate from the chopper settings

    # Fit transmission over the wavelength band corresponding to the lead pulse. The output of the fit
    # (lead_fit_output) has attributes 'fitted_workspace' and 'mantid_fit_output'
    lead_fit_output = fit_band(input_workspace, wavelength_bands.lead, fit_function=fit_function,
                               output_workspace=unique_workspace_dundername())  # band from lead pulse

    # Fit transmission over the wavelength band corresponding to the skipped pulse, if running in skip-frame mode
    if wavelength_bands.skip is not None:
        # skip_fit_output is a namedtuple with attributes 'fitted_workspace' and 'mantid_fit_output'
        skip_fit_output = fit_band(input_workspace, wavelength_bands.skip, fit_function=fit_function,
                                   output_workspace=unique_workspace_dundername())  # skipped pulse

    # The overall fitted transmission workspace is either the transmission over the wavelength range of the lead
    # pulse or the sum of the fitted transmissions over the lead and skip pulses
    if wavelength_bands.skip is None:
        fitted_workspace = CloneWorkspace(lead_fit_output.fitted_workspace, OutputWorkspace=output_workspace)
    else:
        fitted_workspace = Plus(LHSWorkspace=lead_fit_output.fitted_workspace,
                                RHSWorkspace=skip_fit_output.fitted_workspace,
                                OutputWorkspace=output_workspace)

    # Return a dictionary containing comprehensive output from the fit(s)
    r = dict(transmission=fitted_workspace,  # overall fitted transmission workspace
             lead_transmission=lead_fit_output.fitted_workspace,  # fitted transmission workspace over the lead pulse
             lead_mantid_fit=lead_fit_output.mantid_fit_output)  # comprehensive results from Mantid's Fit algorithm

    if wavelength_bands.skip is None:
        r.update(dict(skip_transmission=None, skip_mantid_fit=None))
    else:
        r.update(dict(skip_transmission=skip_fit_output.fitted_workspace,
                      skip_mantid_fit=skip_fit_output.mantid_fit_output))

    return r


def insert_fitted_values(input_workspace, mantid_fit_output,
                         lower_wavelength_boundary, upper_wavelength_boundary):
    r"""
    Substitute transmission raw values with fitted transmission values over
    a range of wavelengths, and set zero elsewhere.

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Workspace to subtitute raw with fitted transmission values
    mantid_fit_output: namedtuple
        Return value after running Mantid's Fit algorithm
    lower_wavelength_boundary: float
        Lower wavelength boundary for the range of fitted transmission values
    upper_wavelength_boundary: float
        Upper wavelength boundary for the range of fitted transmission values
    """
    fitted = mtd[str(input_workspace)]
    y = fitted.dataY(0)
    fitted_values = np.zeros(len(y))
    idx = list(range(lower_wavelength_boundary, upper_wavelength_boundary))
    fitted_values[idx] = mantid_fit_output.OutputWorkspace.dataY(1)
    fitted.dataY(0)[:] = fitted_values


def insert_fitted_errors(input_workspace, mantid_fit_output,
                         lower_wavelength_boundary, upper_wavelength_boundary):
    r"""
    Substitute raw errors with errors derived from the model transmission
    over a wavelength band. Substitute with zero errors elsewhere

    Errors are calculated using the errors in the fitting parameters of the transmission model.
    For instance, the errors in the slope and intercept of a linear model.
    We employ numerical derivative of the fit function with respect to the fitting parameters

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Workspace to contain the fitted error transmission values
    mantid_fit_output: namedtuple
        Return value of Mantid's Fit algorithm
    lower_wavelength_boundary: float
        Lower wavelength boundary of the range of fitted transmission values
    upper_wavelength_boundary: float
        Upper wavelength boundary of the range of fitted transmission values
    """
    fit_function = mantid_fit_output.Function
    wavelengths = mantid_fit_output.OutputWorkspace.dataX(0)
    if len(wavelengths) == len(mantid_fit_output.OutputWorkspace.dataY(0)) + 1:
        wavelengths = (wavelengths[: -1] + wavelengths[1:]) / 2  # dealing with histogram data

    # Iterate over the fitting parameters, calculating their numerical derivatives and their error contributions
    fit_errors = np.zeros(len(wavelengths))
    parameter_table = mantid_fit_output.OutputParameters
    for row_index in range(parameter_table.rowCount() - 1):  # last row is the Chi-square, thus we exclude it
        row = parameter_table.row(row_index)
        parameter_name, parameter_value, parameter_error = [row[key] for key in ('Name', 'Value', 'Error')]
        # evaluate the fit function by increasing the parameter value
        fit_function[parameter_name] += parameter_error  # slightly change the parameter's value
        evaluation_plus = fit_function(wavelengths)  # evaluate function at the domain
        # evaluate the fit function by decreasing the parameter value
        fit_function[parameter_name] -= 2 * parameter_error
        evaluation_minus = fit_function(wavelengths)
        numerical_derivative = (evaluation_plus - evaluation_minus) / (2 * parameter_error)
        fit_errors += (numerical_derivative * parameter_error)**2  # error contribution from this parameter
        fit_function[parameter_name] += parameter_error  # restore the original value
    fit_errors = np.sqrt(fit_errors)

    # Insert errors
    fitted_errors = np.zeros(len(mtd[input_workspace].dataE(0)))
    idx = list(range(lower_wavelength_boundary, upper_wavelength_boundary))
    fitted_errors[idx] = fit_errors
    mtd[input_workspace].dataE(0)[:] = fitted_errors
