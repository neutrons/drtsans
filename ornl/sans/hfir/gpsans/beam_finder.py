from mantid.simpleapi import (MaskBTP, FindCenterOfMassPosition)
from mantid.kernel import logger


def direct_beam_center(input_ws, tubes_to_mask=None):
    '''
    Return beam center x, y in meters
    '''
    if tubes_to_mask is not None:
        MaskBTP(Workspace=input_ws, Tube=tubes_to_mask)

    center = FindCenterOfMassPosition(InputWorkspace=input_ws.name())
    center_x, center_y = center
    logger.information("Found beam position: X={:.3} m, Y={:.3} m.".format(
        center_x, center_y))
    return center_x, center_y
