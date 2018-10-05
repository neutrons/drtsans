
from mantid.simpleapi import (
    LoadEventNexus, SANSMaskDTP, Integration, FindCenterOfMassPosition)


def direct_beam_center(filename, tubes_to_mask=None):
    ws = LoadEventNexus(filename)
    if tubes_to_mask is not None:
        SANSMaskDTP(InputWorkspace=ws, Tube=tubes_to_mask)
    # Flatten TOF
    ws_flattened = Integration(InputWorkspace=ws)
    center = FindCenterOfMassPosition(InputWorkspace=ws_flattened)
    center_x, center_y = center
    return center_x, center_y
