from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import (Integration, FindCenterOfMassPosition)


def direct_beam_center(input_ws, method='center_of_mass'):
    r"""
    Calculate coordinates of beam impinging on the detector

    Parameters
    ----------
    input_ws: EventsWorkspace
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
    ws_flattened = Integration(InputWorkspace=input_ws)
    return method_to_alg[method](InputWorkspace=ws_flattened)
