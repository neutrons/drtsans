from __future__ import print_function

from mantid.simpleapi import (mtd, LoadAscii, ConvertToHistogram,
                              RebinToWorkspace, NormaliseToUnity, Divide,
                              NormaliseByCurrent, Multiply, DeleteWorkspace,
                              ConvertToDistribution)

from ornl.settings import (unique_workspace_dundername as uwd)


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


def load_beam_flux_file(flux, ws_reference=None, output_workspace=None):
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
    if output_workspace is None:
        output_workspace = uwd()  # make a hidden workspace

    LoadAscii(Filename=flux, Separator="Tab", Unit="Wavelength",
              OutputWorkspace=output_workspace)
    ConvertToHistogram(InputWorkspace=output_workspace,
                       OutputWorkspace=output_workspace)
    ConvertToDistribution(Workspace=output_workspace)
    NormaliseToUnity(InputWorkspace=output_workspace,
                     OutputWorkspace=output_workspace)
    if ws_reference is not None:
        RebinToWorkspace(WorkspaceToRebin=output_workspace,
                         WorkspaceToMatch=ws_reference,
                         OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def normalise_by_proton_charge_and_flux(input_workspace, flux,
                                        output_workspace=None):
    r"""Normalises ws to proton and measured flux

    Parameters
    ----------
    input_workspace : MatrixWorkspace
        Workspace to be normalised, rebinned in wavelength.
    flux : Workspace
        Measured beam flux file ws, usually the output of `load_beam_flux_file`

    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    # Normalise by the flux
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=flux,
           OutputWorkspace=output_workspace)
    # Normalize by Proton charge
    NormaliseByCurrent(InputWorkspace=output_workspace,
                       OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def normalise_by_flux(ws, flux, output_workspace=None):
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
    if output_workspace is None:
        output_workspace = uwd()  # output to hidden workspace
    w_flux = load_beam_flux_file(flux, ws_reference=ws)
    normalise_by_proton_charge_and_flux(ws, w_flux,
                                        output_workspace=output_workspace)
    w_flux.delete()

    return mtd[output_workspace]
