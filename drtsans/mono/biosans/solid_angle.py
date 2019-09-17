from drtsans import solid_angle_correction


def apply_solid_angle_correction_main_detector(input_workspace):
    """ Apply solid angle correction for the main detector

    Parameters
    ----------
    input_workspace : MatrixWorkspace, str
        The input workspace name or itself

    Returns
    -------
    MatrixWorkspace
        The input workspace corrected for solid angle
    """

    return solid_angle_correction(
        input_workspace, detector_type='VerticalTube')


def apply_solid_angle_correction_wing_detector(input_workspace):
    """ Apply solid angle correction for the wing detector

    Parameters
    ----------
    input_workspace : MatrixWorkspace, str
        The input workspace name or itself

    Returns
    -------
    MatrixWorkspace
        The input workspace corrected for solid angle
    """

    return solid_angle_correction(
        input_workspace, detector_type='VerticalWing')
