import drtsans.beam_finder as bf

__all__ = ['center_detector', 'find_beam_center']


def find_beam_center(input_workspace, method='center_of_mass', mask=None, **kwargs):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    based on drtsans.beam_finder.find_beam_center

    Parameters
    ----------
    input_workspace: str, ~mantid.api.Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str, ``MaskWorkspace``
        Path to mask file, or ``MaskWorkspace`` object
    kwargs: dict
        Parameters to be passed to the method to calculate the center or to MaskBTP or 'panel'
        'panel' is either 'front' or 'back' to mask a whole panel

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center, always in meters.
    """
    return bf.find_beam_center(input_workspace, method, mask, **kwargs)


def center_detector(input_workspace, center_x, center_y):
    """Center the detector

    based on drtsans.beam_finder.center_detector

    Parameters
    ----------
    input_workspace : Workspace2D, str
        The workspace to be centered
    center_x : float
        in meters
    center_y : float
        in meters
    """
    bf.center_detector(input_workspace, center_x, center_y, component='detector1')
