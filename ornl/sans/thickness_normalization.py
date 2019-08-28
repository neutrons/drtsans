
def normalize_by_thickness(ws, thickness=None):
    r"""Normalize input workspace by thickness

    Parameters
    ----------
    ws: Input workspace
    thickness: thickness of the sample in centimeters.
        If None, should obtain thickness from log values in the workspace

    Returns
    -------
    normalized workspace
    """
    if thickness is None:
        raise NotImplementedError
    from mantid.simpleapi import Divide
    return ws/thickness
