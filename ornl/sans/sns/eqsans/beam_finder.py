
from mantid.simpleapi import (
    Load, SANSMaskDTP, Integration, FindCenterOfMassPosition)


def direct_beam_center(filename, tubes_to_mask):
    ws = Load(filename)
    SANSMaskDTP(InputWorkspace=ws, Tube=tubes_to_mask)
    # Flatten TOF
    ws_flattened = Integration(InputWorkspace=ws)
    center = FindCenterOfMassPosition(InputWorkspace=ws_flattened)
    center_x, center_y = center
    return center_x, center_y
