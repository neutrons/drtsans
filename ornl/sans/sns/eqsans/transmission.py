from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.simpleapi import (Fit, CloneWorkspace, RenameWorkspace)

from ornl.settings import namedtuplefy
from ornl.sans.sns.eqsans.correct_frame import transmitted_bands
from ornl.sans.transmission import _calculate_radius_from_input_ws


def beam_radius(ws):
    r"""
    Calculate the beam radius impinging on the detector bank.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace, contains all necessary info in the logs

    Returns
    -------
    float
        Estimated beam radius
    """
    logs = dict(sample_aperture_diameter_log='sample-aperture-diameter',
                source_aperture_diameter_log='source-aperture-diameter',
                sdd_log='sample-detector-distance',
                ssd_log='source-aperture-sample-distance')
    return _calculate_radius_from_input_ws(ws, **logs)


def insert_fitted_values(mfit, fitted, low_b, up_b):
    r"""
    Substituted raw with fitted transmission values

    Parameters
    ----------
    mfit: namedtuple
        Return value of Mantid's Fit algorithm
    fitted: MatrixWorkspace
        Workspace to contain the fitted transmission values
    low_b: float
        Lower wavelength boundary of the range of fitted transmission values
    up_b: float
        Upper wavelength boundary of the range of fitted transmission values
    """
    y = fitted.dataY(0)
    ins = np.zeros(len(y))
    idx = list(range(low_b, up_b))
    ins[idx] = mfit.OutputWorkspace.dataY(1)
    fitted.dataY(0)[:] = ins


def insert_fitted_errors(mfit, fitted, low_b, up_b):
    r"""
    Substitute raw errors with errors derived from the model transmission.

    Errors are calculated using the errors in the fitting parameters of the
    transmission model. For instance, the errors in the slope and intercept
    of a linear model.

    Parameters
    ----------
    mfit: namedtuple
        Return value of Mantid's Fit algorithm
    fitted: MatrixWorkspace
        Workspace to contain the fitted error transmission values
    low_b: float
        Lower wavelength boundary of the range of fitted transmission values
    up_b: float
        Upper wavelength boundary of the range of fitted transmission values
    """
    # Estimate errors using the numerical derivative of the function with
    # respect to the fitting parameters
    f = mfit.Function
    x = mfit.OutputWorkspace.dataX(0)
    if len(x) == len(mfit.OutputWorkspace.dataY(0)) + 1:
        x = (x[: -1] + x[1:]) / 2  # dealing with histogram data
    e = np.zeros(len(x))
    p_table = mfit.OutputParameters
    for i in range(p_table.rowCount() - 1):
        row = p_table.row(i)
        p_n, p_e = row['Name'], row['Error']
        f[p_n] = f[p_n] + p_e  # slightly change the parameter's value
        d = f(x)  # evaluate function at the domain
        f[p_n] = f[p_n] - 2 * p_e
        d = (d - f(x)) / (2 * p_e)  # numerical derivative with respect to p_n
        e += np.abs(d) * p_e**2  # error contribution
        f[p_n] += p_e
    e = np.sqrt(e)

    # Insert errors
    ins = np.zeros(len(fitted.dataE(0)))
    idx = list(range(low_b, up_b))
    ins[idx] = e
    fitted.dataE(0)[:] = ins


@namedtuplefy
def fit_band(raw, band, func, suffix=None):
    r"""
    Fit the wavelength dependence of the raw zero-angle transmission
    values with a function within a wavelength band.

    Parameters
    ----------
    raw: MatrixWorkspace
        Workspace containing the raw transmission
    band: Wband
        Wavelength band over which to carry out the fit
    func: str
        String representation of the fit function. See Mantid's
        `UserFunction` or any of Mantid's fit functions
    suffix: str
        suffix for names of output workspaces.

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - fitted: MatrixWorkspace, transmission values within the band,
            zero elsewhere
        - mfit: namedtuple, output when calling Mantid's Fit algorithm
    """
    # Carry out the fit only on the wavelength band
    x = raw.dataX(0)
    y = raw.dataY(0)
    bi = np.where((x >= band.min) & (x < band.max))[0]
    min_y = 1e-3 * np.mean(y[bi[:-1]])  # 1e-3 pure heuristics

    # Find wavelength range with non-zero intensities. Care with boundaries
    i = 0
    while x[i] < band.min or y[i] < min_y:
        i += 1
    lower_bin_boundary = i
    while i < len(y) and x[i] < band.max and y[i] > min_y:
        i += 1
    upper_bin_boundary = i
    start_x, end_x = x[lower_bin_boundary], x[upper_bin_boundary - 1]
    mfit = Fit(Function=func, InputWorkspace=raw.name(), WorkspaceIndex=0,
               StartX=start_x, EndX=end_x, Output=raw.name() + '_fit')

    # Insert the fitted band into the wavelength range of raw
    name = '{}_fitted_{}'.format(raw.name(), suffix)
    fitted = CloneWorkspace(raw, OutputWorkspace=name)
    insert_fitted_values(mfit, fitted, lower_bin_boundary, upper_bin_boundary)
    insert_fitted_errors(mfit, fitted, lower_bin_boundary, upper_bin_boundary)
    return dict(fitted=fitted, mfit=mfit)


@namedtuplefy
def fit_raw(raw, fitted, func='name=UserFunction,Formula=a*x+b'):
    r"""
    Fit the wavelength dependence of the raw zero-angle transmission
    values with a model.

    If working in frame skipping mode, apply the fit separately to the
    wavelength bands of the lead and skipped pulses.

    Parameters
    ----------
    raw: MatrixWorkspace
        Workspace containing the raw transmission
    fitted: str
        Name of the output workspace containing the fitted transmission
        values
    func: str
        String representation of the fit function. See Mantid's
        `UserFunction` or any of Mantid's fit functions

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - fit: workspace containing the fitted transmission values and errors
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

    # Fit only over the range of the transmitted wavelength band(s)
    bands = transmitted_bands(raw)
    fit_lead = fit_band(raw, bands.lead, func, 'lead')  # band from lead pulse
    fitted_ws = fit_lead.fitted
    if bands.skip is not None:
        fit_skip = fit_band(raw, bands.skip, func, 'skip')  # skipped pulse
        fitted_ws += fit_skip.fitted
    fitted_ws = RenameWorkspace(fitted_ws, OutputWorkspace=fitted,
                                RenameMonitors=False)
    # dictionary to return
    r = dict(transmission=fitted_ws,
             lead_fit=fit_lead.fitted, lead_mfit=fit_lead.mfit)
    if bands.skip is not None:
        r.update(dict(skip_fit=fit_skip.fitted, skip_mfit=fit_skip.mfit))
    else:
        r.update(dict(skip_fit=None, skip_mfit=None))
    return r
