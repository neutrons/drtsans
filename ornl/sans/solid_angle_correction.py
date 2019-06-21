from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from mantid.simpleapi import SolidAngle, ReplaceSpecialValues


def solid_angle_correction(input_workspace, detector_type):
    solid_angle_ws = SolidAngle(InputWorkspace=input_workspace,
                                    Method=detector_type)
    output_workspace = input_workspace / solid_angle_ws
    ReplaceSpecialValues(InputWorkspace=output_workspace,
                         OutputWorkspace=output_workspace, NaNValue=0.,
                         InfinityValue=0.)
    return output_workspace
