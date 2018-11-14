
from mantid.simpleapi import (
    SANSMaskDTP, Integration, FindCenterOfMassPosition)


def direct_beam_center(input_ws, tubes_to_mask=None):
    if tubes_to_mask is not None:
        SANSMaskDTP(InputWorkspace=input_ws, Tube=tubes_to_mask)
    # Flatten TOF
    ws_flattened = Integration(InputWorkspace=input_ws)
    center = FindCenterOfMassPosition(InputWorkspace=ws_flattened)
    center_x, center_y = center
    return center_x, center_y
