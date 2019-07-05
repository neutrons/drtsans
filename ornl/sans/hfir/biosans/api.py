from ornl.sans import solid_angle_correction


def apply_solid_angle_correction(ws):
    """
        Apply solid angle correction
    """
    return solid_angle_correction(ws, detector_type='VerticalTube')
