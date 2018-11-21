
from mantid.simpleapi import (
    LoadAscii, ConvertToHistogram, RebinToWorkspace, NormaliseToUnity, Divide,
    NormaliseByCurrent, Multiply, DeleteWorkspace, Integration)


def load_beam_flux_file(file_path, ws_reference=None):
    '''Loads the ascii beam flux file

    Parameters
    ----------
    ws_reference : Workspace
        The reference workspace to rebin the flux to. If none, does not
        rebin the data.
    '''

    ws = LoadAscii(Filename=file_path, Separator="Tab", Unit="Wavelength")
    ws = ConvertToHistogram(InputWorkspace=ws)
    if ws_reference is not None:
        ws = RebinToWorkspace(WorkspaceToRebin=ws,
                              WorkspaceToMatch=ws_reference)
    ws = NormaliseToUnity(InputWorkspace=ws)
    return ws


def time(ws_input, t_frame, t_low_cut, t_high_cut, bin_width, wavelength_min, wavelength_max):
    '''Time normalization. It is only used for dark current.
    It does not make sense to use it for other files
    
    
    Parameters
    ----------
    ws_input : Workspace
        The dark current ws
    t_frame : float
        
    t_low_cut : [type]
        
    t_high_cut : [type]
        [description]
    bin_width : [type]
        [description]
    wavelength_min : [type]
        [description]
    wavelength_max : [type]
        [description]
    
    '''

    # Not done!!!
    


    # First remove the time component
    # If it's already in this format the ingration does nothing
    ws_dark_current = Integration(InputWorkspace=ws_dark_current)

    #TODO
    # We think that uses information from a sample (or other dataset used to 
    # subtract the DC later) 





def monitor(ws_input, ws_monitor, ws_flux_to_monitor_ratio):
    '''Monitor normalisation
    
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
    '''

    ws_tmp = Multiply(LHSWorkspace=ws_monitor,
                      RHSWorkspace=ws_flux_to_monitor_ratio)
    ws = Divide(LHSWorkspace=ws_input, RHSWorkspace=ws_tmp)
    DeleteWorkspace(ws_tmp)
    return ws


def proton_charge_and_flux(ws_input, ws_beam_flux):
    """Normalises ws to proton and measured flux

    Parameters
    ----------
    ws_input : Workspace
        The workspace to be normalised

    ws_beam_flux : Workspace
        The measured beam flux file ws

    """
    # Normalise by the flux
    ws = Divide(LHSWorkspace=ws_input, RHSWorkspace=ws_beam_flux)
    # Proton charge
    ws = NormaliseByCurrent(InputWorkspace=ws)
    return ws
