from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from mantid.simpleapi import SolidAngle, ReplaceSpecialValues
from ornl.settings import optional_output_workspace


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
                                Method=detector_type)
    output_workspace = input_workspace / solid_angle_ws
    ReplaceSpecialValues(InputWorkspace=output_workspace,
                         OutputWorkspace=output_workspace, NaNValue=0.,
                         InfinityValue=0.)
    return output_workspace
