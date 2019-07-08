from __future__ import print_function

from mantid.simpleapi import (
    LoadAscii, ConvertToHistogram, RebinToWorkspace, NormaliseToUnity, Divide,
    NormaliseByCurrent, Multiply, DeleteWorkspace, ConvertToDistribution)

from ornl.settings import (optional_output_workspace,
                           unique_workspace_dundername as uwd)


__all__ = ['normalise_by_flux', ]


def monitor(ws_input, ws_monitor, ws_flux_to_monitor_ratio):
    r"""Monitor normalisation

    Parameters
    ----------
    ws_input : Workspace
        The workspace to be normalised
    ws_monitor : Workspace
        The workspace with the monitor count
    ws_flux_to_monitor_ratio : Workspace
        Pre-mesured flux-to-monitor ratio spectrum

    Returns
    -------
    Workspace
        [description]
    """

    ws_tmp = Multiply(LHSWorkspace=ws_monitor,
                      RHSWorkspace=ws_flux_to_monitor_ratio)
    ws = Divide(LHSWorkspace=ws_input, RHSWorkspace=ws_tmp)
    DeleteWorkspace(ws_tmp)
    return ws


@optional_output_workspace
def load_beam_flux_file(flux, ws_reference=None):
    r"""

    Loads the beam flux file and converts to a wavelength
    normalized probability distribution.

    Parameters
    ----------
    flux: str
        Path to file with the wavelength distribution of the neutron
        flux. Loader is Mantid `LoadAscii` algorithm.

    ws_reference : Workspace
        Workspace to rebin the flux to. If None, no rebin is performed
    """

    ws = LoadAscii(Filename=flux, Separator="Tab", Unit="Wavelength",
                   OutputWorkspace=uwd())
    ws = ConvertToHistogram(InputWorkspace=ws,
                            OutputWorkspace=ws.name())
    ConvertToDistribution(ws)
    ws = NormaliseToUnity(ws, OutputWorkspace=ws.name())
    if ws_reference is not None:
        ws = RebinToWorkspace(WorkspaceToRebin=ws,
                              WorkspaceToMatch=ws_reference,
                              OutputWorkspace=ws.name())
    return ws


@optional_output_workspace
def normalise_by_proton_charge_and_flux(ws, flux):
    r"""Normalises ws to proton and measured flux

    Parameters
    ----------
    ws : MatrixWorkspace
        Workspace to be normalised, rebinned in wavelength.
    flux : Workspace
        Measured beam flux file ws, usually the output of `load_beam_flux_file`

    """
    # Normalise by the flux
    w = Divide(LHSWorkspace=ws, RHSWorkspace=flux, OutputWorkspace=uwd())
    # Normalize by Proton charge
    w = NormaliseByCurrent(w, OutputWorkspace=w.name())
    return w


@optional_output_workspace
def normalise_by_flux(ws, flux):
    r"""
    Normalize counts by flux wavelength distribution and proton charge.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace, binned in wavelength
    flux: str
        path to file containing the wavelength distribution
        of the neutron flux.

    Returns
    -------
    MatrixWorkspace
    """
    w_flux = load_beam_flux_file(flux, ws_reference=ws, output_workspace=uwd())
    return normalise_by_proton_charge_and_flux(ws, w_flux,
                                               output_workspace=uwd())
