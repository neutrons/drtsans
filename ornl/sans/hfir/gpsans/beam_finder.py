
from mantid.simpleapi import (
    LoadSpice2D, SANSMaskDTP, FindCenterOfMassPosition)


def direct_beam_center(filename, tubes_to_mask=None):
    '''
    Return beam center x, y in meters
    '''
    ws = LoadSpice2D(filename)
    if tubes_to_mask is not None:
        SANSMaskDTP(InputWorkspace=ws.OutputWorkspace, Tube=tubes_to_mask)

    center = FindCenterOfMassPosition(InputWorkspace=ws.OutputWorkspace)
    center_x, center_y = center

    return center_x, center_y

