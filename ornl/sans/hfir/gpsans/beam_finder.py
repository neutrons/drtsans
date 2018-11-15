from mantid.simpleapi import (
    SANSMaskDTP, FindCenterOfMassPosition)


def direct_beam_center(input_ws, tubes_to_mask=None):
    '''
    Return beam center x, y in meters
    '''
    if tubes_to_mask is not None:
        SANSMaskDTP(InputWorkspace=input_ws, Tube=tubes_to_mask)

    center = FindCenterOfMassPosition(InputWorkspace=input_ws)
    center_x, center_y = center

    return center_x, center_y
