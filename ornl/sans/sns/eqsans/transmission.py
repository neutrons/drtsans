from __future__ import (absolute_import, division, print_function)

import numpy as np
from ornl.settings import namedtuplefy
from ornl.sans.sns.eqsans.correct_frame import transmitted_bands
from mantid.simpleapi import Fit, CloneWorkspace, RenameWorkspace


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
    bi = np.where((x >= band.min) & (x <= band.max))
    mfit = Fit(func, raw, WorkspaceIndex=0, StartX=bi[0], EndX=bi[-1],
               Output=raw.name + '_fit')

    # Insert the fitted band into the wavelength range of raw
    name = '{}_fitted_{}'.format(raw.name(), suffix)
    fitted = CloneWorkspace(raw, OutputWorkspace=name)
    ins = np.zeros(len(y))
    shift = 0 if len(x) == len(y) else 1  # dealing with histograms
    ins[bi + shift] = fitted.OutputWorkspace.dataY(1)
    fitted.dataY(0)[:] = ins

    return dict(fitted=fitted, mfit=mfit)


@namedtuplefy
def fit_raw(raw, fitted, func='name=UserFunction,Formula=a*x+b'):
    r"""
    Fit the wavelength dependence of the raw zero-angle transmission
    values with a function.
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
        - fit: workspace containing only the fitted transmission values
            (errors are those of the raw workspace)
    """

    # Fit only over the range of the transmitted wavelength band(s)
    bands = transmitted_bands(raw)
    fit_lead = fit_band(raw, bands.lead, func, 'lead')  # band from lead pulse
    fitted_ws = fit_lead.fitted
    if bands.skip is not None:
        fit_skip = fit_band(raw, bands.skip, func, 'skip')  # skipped pulse
        fitted_ws += fit_skip.fitted
    fitted_ws = RenameWorkspace(fitted_ws, fitted)

    return dict(fit=fitted_ws)
