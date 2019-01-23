from __future__ import (absolute_import, division, print_function)


from mantid.simpleapi import (Integration, FindCenterOfMassPosition)
from ornl.settings import unique_workspace_name


def direct_beam_center(input_ws,
                       finder=FindCenterOfMassPosition, finder_kwargs={}):
    r"""
    Calculate coordinates of beam impinging on the detector for a
    direct-beam run

    Parameters
    ----------
    input_ws: EventsWorkspace
        Workspace for the direct beam run (no sample and no sample holder)
    finder: function
        Method to calculate the beam center
    finder_kwargs: dict
        Additional options for the finder method as a python dictionary

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (in meters)
    """
    ws_flattened = Integration(InputWorkspace=input_ws,
                               OutputWorkspace=unique_workspace_name())
    return finder(InputWorkspace=ws_flattened, **finder_kwargs)
