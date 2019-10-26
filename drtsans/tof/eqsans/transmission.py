import numpy as np
from mantid.api import mtd
from mantid.simpleapi import (Fit, CloneWorkspace, Plus, SaveNexus)

from drtsans.settings import (namedtuplefy, unique_workspace_dundername as uwd)
from drtsans.transmission import calculate_transmission as calculate_raw_transmission
from drtsans.tof.eqsans.correct_frame import transmitted_bands
from drtsans.tof.eqsans.geometry import sample_aperture_diameter, source_aperture_diameter
from drtsans.geometry import sample_detector_distance, source_sample_distance


# Symbols to be exported to the eqsans namespace
__all__ = ['beam_radius', 'calculate_transmission']


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
        output_workspace = uwd()
    transmission_values = calculate_raw_transmission(input_sample, input_reference, radius=radius,
                                                     radius_unit=radius_unit, output_workspace=output_workspace)
    if bool(fit_function) is True:
        transmission_values = fit_raw(transmission_values, func=fit_function).transmission
    return transmission_values


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
def fit_band(input_workspace, band, func='name=UserFunction,Formula=a*x+b',
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
    func: str
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
        - fitted: MatrixWorkspace, transmission values within the band,
            zero elsewhere. The name of this workspace is `output_workspace`
        - mfit: namedtuple, output of calling Mantid's Fit algorithm
    """
    if output_workspace is None:
        output_workspace = uwd()
    ws = mtd[str(input_workspace)]
    # Carry out the fit only over the wavelength band with sensible intensities
    x = ws.dataX(0)
    y = ws.dataY(0)
    band_indexes = np.where((x >= band.min) & (x < band.max))[0]
    min_y = 1e-3 * np.mean(y[band_indexes[:-1]])  # 1e-3 pure heuristics

    # Find wavelength range with non-zero intensities. Care with boundaries
    # that have intensities largely deviating from the expected intensities
    i = 0
    while x[i] < band.min or y[i] < min_y:
        i += 1
    lower_bin_boundary = i
    while i < len(y) and x[i] < band.max and y[i] > min_y:
        i += 1
    upper_bin_boundary = i
    start_x, end_x = x[lower_bin_boundary], x[upper_bin_boundary - 1]
    SaveNexus(ws, '/tmp/junk.nxs')
    mfit = Fit(Function=func, InputWorkspace=ws.name(), WorkspaceIndex=0,
               StartX=start_x, EndX=end_x, Output=uwd())

    # Insert the fitted band into the wavelength range of ws
    fitted = CloneWorkspace(ws, OutputWorkspace=output_workspace)
    insert_fitted_values(output_workspace, mfit,
                         lower_bin_boundary, upper_bin_boundary)
    insert_fitted_errors(output_workspace, mfit,
                         lower_bin_boundary, upper_bin_boundary)
    return dict(fitted=fitted, mfit=mfit)


@namedtuplefy
def fit_raw(input_workspace, func='name=UserFunction,Formula=a*x+b',
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
    func: str
        String representation of the fit function. See Mantid's
        `UserFunction` or any of Mantid's fit functions

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - transmission: workspace containing the fitted transmission
               values and errors.
               The name of this workspace is the value of `output_workspace`
        - lead_fit: workspace containing the fitted transmission values and
            errors of the lead pulse
        - lead_mfit: return value of running Mantid's Fit algorithm when
            fitting the raw transmission over the lead pulse wavelength range
        - skip_fit: workspace containing the fitted transmission values and
            errors of the skip pulse. None if not working in frame skipping
            mode
        - skip_mfit: return value of running Mantid's Fit algorithm when
            fitting the raw transmission over the skip pulse wavelength range
            None f not working in frame skipping mode
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    raw = mtd[str(input_workspace)]
    #
    # Fit only over the range of the transmitted wavelength band(s)
    #
    bands = transmitted_bands(raw)
    #
    # Fit the lead pulse
    #
    fit_lead = fit_band(raw, bands.lead, func=func,
                        output_workspace=uwd())  # band from lead pulse
    #
    # Fit the skipped pulse, if running in skip-frame mode
    #
    if bands.skip is not None:
        fit_skip = fit_band(raw, bands.skip, func=func,
                            output_workspace=uwd())  # skipped pulse
        fitted_ws = Plus(LHSWorkspace=fit_lead.fitted,
                         RHSWorkspace=fit_skip.fitted,
                         OutputWorkspace=output_workspace)
    else:
        fitted_ws = CloneWorkspace(fit_lead.fitted,
                                   OutputWorkspace=output_workspace)
    # dictionary to return
    r = dict(transmission=fitted_ws,
             lead_fit=fit_lead.fitted, lead_mfit=fit_lead.mfit)
    if bands.skip is not None:
        r.update(dict(skip_fit=fit_skip.fitted, skip_mfit=fit_skip.mfit))
    else:
        r.update(dict(skip_fit=None, skip_mfit=None))
    return r


def insert_fitted_values(input_workspace, mantid_fit_output,
                         lower_wavelength_boundary, upper_wavelength_boundary):
    r"""
    Substitute transmission raw values with fitted transmission values over
    a range of vawelengths, and set zero elsewhere.

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

    Errors are calculated using the errors in the fitting parameters of the
    transmission model. For instance, the errors in the slope and intercept
    of a linear model.

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
    fitted = mtd[str(input_workspace)]
    #
    # Estimate errors using the numerical derivative of the fit function with
    # respect to the fitting parameters
    #
    f = mantid_fit_output.Function
    x = mantid_fit_output.OutputWorkspace.dataX(0)
    if len(x) == len(mantid_fit_output.OutputWorkspace.dataY(0)) + 1:
        x = (x[: -1] + x[1:]) / 2  # dealing with histogram data
    e = np.zeros(len(x))
    p_table = mantid_fit_output.OutputParameters
    # Iterate over the fitting parameters,calculating the numerical derivative
    for i in range(p_table.rowCount() - 1):
        row = p_table.row(i)
        p_n, p_e = row['Name'], row['Error'] / 2.0
        f[p_n] = f[p_n] + p_e  # slightly change the parameter's value
        d = f(x)  # evaluate function at the domain
        f[p_n] = f[p_n] - 2 * p_e
        d = (d - f(x)) / (2 * p_e)  # numerical derivative with respect to p_n
        e += np.abs(d) * p_e**2  # error contribution from this parameter
        f[p_n] += p_e  # restore the original value
    e = np.sqrt(e)

    # Insert errors
    fitted_errors = np.zeros(len(fitted.dataE(0)))
    idx = list(range(lower_wavelength_boundary, upper_wavelength_boundary))
    fitted_errors[idx] = e
    fitted.dataE(0)[:] = fitted_errors
