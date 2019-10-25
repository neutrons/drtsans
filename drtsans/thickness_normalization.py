
def normalize_by_thickness(ws, thickness):
    r"""Normalize input workspace by thickness

    Parameters
    ----------
    ws: Input workspace
    thickness: thickness of the sample in centimeters.

    Returns
    -------
    normalized workspace
    """
    if thickness <= 0.0:
        msg = 'Sample thickness should be positive. Got {}'.format(thickness)
        raise ValueError(msg)
    return ws/thickness
