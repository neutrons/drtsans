from drtsans import solid_angle_correction

__all__ = ['solid_angle_correction_main_detector', 'solid_angle_correction_wing_detector']


def solid_angle_correction_main_detector(input_workspace):
    """ Apply :func:`drtsans.solid_angle_correction` to the main detector

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace, str
        The input workspace name or itself

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The input workspace corrected for solid angle
    """

    return solid_angle_correction(input_workspace, detector_type='VerticalTube')


def solid_angle_correction_wing_detector(input_workspace):
    """ Apply :func:`drtsans.solid_angle_correction` to the wing detector

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace, str
        The input workspace name or itself

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The input workspace corrected for solid angle
    """

    return solid_angle_correction(
        input_workspace, detector_type='VerticalWing')
