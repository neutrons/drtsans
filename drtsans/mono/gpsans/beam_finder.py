from mantid.simpleapi import (MaskBTP, FindCenterOfMassPosition)
from mantid.kernel import logger


def find_beam_center(input_ws, method='center_of_mass', mask=None, **kwargs):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    Parameters
    ----------
    input_workspace: str, Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str or list
        Tubes to be masked
    kwargs: dict
        Parameters to be passed to the method to calculate the center

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters)
    """
    if method != 'center_of_mass':
        raise NotImplementedError(f'{method} is not implemented')
    # TODO: use apply_mask instead
    if mask is not None:
        MaskBTP(Workspace=input_ws, Tube=mask)

    center = FindCenterOfMassPosition(InputWorkspace=input_ws.name(), **kwargs)
    center_x, center_y = center
    logger.information("Found beam position: X={:.3} m, Y={:.3} m.".format(
        center_x, center_y))
    return center_x, center_y
