from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from mantid.simpleapi import (Divide, SolidAngle,
                              ReplaceSpecialValues)
from ornl.settings import (optional_output_workspace,
                           unique_workspace_dundername as uwd)


@optional_output_workspace
def solid_angle_correction(input_workspace, detector_type='Rectangle'):
    r"""
    The algorithm calculates solid angles from the sample position of
    the input workspace for all of the spectra selected. The output workspace
    is the input divided by the solid angle.

    Parameters
    __________

    input_workspace: MatrixWorkspace

    detector_type: Select the method to calculate the Solid Angle. Allowed
    values: [‘GenericShape’, ‘Rectangle’, ‘VerticalTube’, ‘HorizontalTube’,
    ‘VerticalWing’, ‘HorizontalWing’]

    Returns
    _______
    MatrixWorkspace with the solid angle correction applied.

    """
    solid_angle_ws = SolidAngle(InputWorkspace=input_workspace,
                                OutputWorkspace=uwd(), 
                                Method=detector_type)
    output_workspace = Divide(LHSWorkspace=input_workspace,
                              RHSWorkspace=solid_angle_ws,
                              OutputWorkspace=uwd())
    ReplaceSpecialValues(InputWorkspace=output_workspace,
                         OutputWorkspace=output_workspace, NaNValue=0.,
                         InfinityValue=0.)
    return output_workspace
