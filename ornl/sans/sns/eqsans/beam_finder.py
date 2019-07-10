from mantid.simpleapi import (Integration, FindCenterOfMassPosition)
from ornl.settings import unique_workspace_dundername as uwd


def direct_beam_center(ws, method='center_of_mass'):
    r"""
    Calculate absolute coordinates of beam impinging on the detector

    Parameters
    ----------
    ws: EventsWorkspace
        Workspace for the direct beam run (no sample and no sample holder)
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters)
    """
    method_to_alg = dict(center_of_mass=FindCenterOfMassPosition)
    ws_flattened = Integration(InputWorkspace=ws, OutputWorkspace=uwd())
    return method_to_alg[method](InputWorkspace=ws_flattened)
