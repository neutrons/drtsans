from mantid.simpleapi import (mtd, LoadAscii, ConvertToHistogram,
                              RebinToWorkspace, NormaliseToUnity, Divide,
                              NormaliseByCurrent, ConvertToDistribution,
                              CloneWorkspace, RemoveSpectra, Multiply, Load,
                              DeleteWorkspace)
from ornl import path
from ornl.settings import (unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs

__all__ = ['normalise_by_flux', ]


def load_beam_flux_file(flux, ws_reference=None, output_workspace=None):
    r"""

    Loads the beam flux file and converts to a wavelength
    normalized probability distribution.

    Parameters
    ----------
    flux: str
        Path to file with the wavelength distribution of the neutron
        flux. Loader is Mantid `LoadAscii` algorithm.
    ws_reference : str, MatrixWorkspace
        Workspace to rebin the flux to. If None, no rebin is performed
    output_workspace: str
        Name of the output workspace. If None, a hidden random name
        will be assigned.

    Returns
    -------
    MatrixWorkspace
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
    r"""
    Normalises the input workspace by proton charge and measured flux

    Parameters
    ----------
    input_workspace : str, MatrixWorkspace
        Workspace to be normalised, rebinned in wavelength.
    flux : Workspace
        Measured beam flux file ws, usually the output of `load_beam_flux_file`
    output_workspace : str
        Name of the normalised workspace. If None, the name of the input
        workspace is chosen (the input workspace is overwritten).

    Returns
    -------
    MatrixWorkspace
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


def load_flux_to_monitor_ratio_file(flux, ws_reference=None,
                                    loader_kwargs=dict(),
                                    output_workspace=None):
    r"""
    Loads the flux-to-monitor ratio

    Parameters
    ----------
    flux: str
        Path to file with the flux-to-monitor ratio data. Loader is
        Mantid `LoadAscii` algorithm.
    loader_kwargs: dict
        optional keyword arguments to Mantid's Load algorithm
    output_workspace: str
        Name of the output workspace. If None, a hidden random name
        will be assigned.

    Returns
    -------
    MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = uwd()  # make a hidden workspace

    # Let Mantid figure out what kind file format is the flux file
    Load(Filename=flux, OutputWorkspace=output_workspace, **loader_kwargs)
    ConvertToHistogram(InputWorkspace=output_workspace,
                       OutputWorkspace=output_workspace)
    if ws_reference is not None:
        # RebinToWorkspace akin to average of the ratios within each bin,
        # except of a scaling factor
        RebinToWorkspace(WorkspaceToRebin=output_workspace,
                         WorkspaceToMatch=ws_reference,
                         OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def normalise_by_monitor(input_workspace, flux_to_monitor, monitor_workspace,
                         output_workspace=None):
    r"""
    Normalises the input workspace by monitor count and flux-to-monitor
    ratio.

    Parameters
    ----------
    input_workspace : str, MatrixWorkspace
        Workspace to be normalised, rebinned in wavelength.
    flux_to_monitor : str, MatrixWorkspace
        Flux to monitor ratio. A file path or a workspace resulting from
        calling `load_flux_to_monitor_ratio_file`.
    monitor_workspace : str, MatrixWorkspace
        Counts from the monitor.
    output_workspace : str
        Name of the normalised workspace. If None, the name of the input
        workspace is chosen (the input workspace is overwritten).

    Returns
    -------
    MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    wi = mtd[str(input_workspace)]

    # Check non-skip mode
    if bool(SampleLogs(wi).is_frame_skipping.value) is True:
        msg = 'Normalisation by monitor not possible in frame-skipping mode'
        raise ValueError(msg)

    # Only the first spectrum of the monitor is required
    monitor = uwd()
    RebinToWorkspace(monitor_workspace, input_workspace,
                     OutputWorkspace=monitor)
    excess_idx = range(1, mtd[monitor].getNumberHistograms())
    RemoveSpectra(monitor, WorkspaceIndices=excess_idx,
                  OutputWorkspace=monitor)

    # Elucidate the nature of the flux to monitor input
    flux = uwd()
    if isinstance(flux_to_monitor, str) and path.exists(flux_to_monitor):
        load_flux_to_monitor_ratio_file(flux_to_monitor,
                                        ws_reference=input_workspace,
                                        output_workspace=flux)
    else:
        CloneWorkspace(flux_to_monitor, OutputWorkspace=flux)
        RebinToWorkspace(flux, input_workspace, OutputWorkspace=flux)

    # Generate a normalized beam flux distribution
    Multiply(flux, monitor, OutputWorkspace=flux)
    ConvertToDistribution(Workspace=flux)
    NormaliseToUnity(InputWorkspace=flux, OutputWorkspace=flux)

    # Normalise our input workspace
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=flux,
           OutputWorkspace=output_workspace)

    # Clean the dust balls
    [DeleteWorkspace(name) for name in (flux, monitor)]
    return mtd[output_workspace]


def normalise_by_flux(input_workspace, flux, method='proton charge',
                      monitor_workspace=None, output_workspace=None):
    r"""
    Normalize counts by several methods to estimate the neutron flux.

    Neutron flux can be estimated either with a monitor or with the proton
    charge.

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Input workspace, binned in wavelength
    flux: str
        If `method` is 'proton charge', flux is the path to the file
        containing the wavelength distribution of the neutron flux. If
        `method` is `monitor`, then flux is the path to the file containing
        a pre-measured flux-to-monitor ratio spectrum
    method: str
        Either 'proton charge' or 'monitor'
    monitor_workspace: str, MatrixWorkspace
        Prepared monitor workspace
    output_workspace : str
        Name of the normalised workspace. If None, the name of the input
        workspace is chosen (the input workspace is overwritten).

    Returns
    -------
    MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    # Select the flux file
    flux_loader = {'proton charge': load_beam_flux_file,
                   'monitor': load_flux_to_monitor_ratio_file}
    w_flux = flux_loader[method](flux, ws_reference=input_workspace)
    # Select the normalisation function
    normaliser = {'proton charge': normalise_by_proton_charge_and_flux,
                  'monitor': normalise_by_monitor}
    # Additional arguments specific to the normaliser
    args = {'proton charge': list(),
            'monitor': [monitor_workspace]}
    normaliser[method](input_workspace, w_flux, *args[method],
                       output_workspace=output_workspace)
    # A bit of cleanup
    w_flux.delete()

    return mtd[output_workspace]
