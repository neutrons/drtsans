from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
# https://docs.mantidproject.org/nightly/algorithms/SolidAngle-v1.html
# https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html
from mantid.simpleapi import (DeleteWorkspace, Divide, mtd, SolidAngle,
                              ReplaceSpecialValues)
from ornl.settings import (unique_workspace_dundername as uwd)


def solid_angle_correction(input_workspace, detector_type='VerticalTube',
                           output_workspace=None):
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
    input_workspace = str(input_workspace)
    solid_angle_ws = uwd()  # temporary workspace
    if output_workspace is None:
        output_workspace = input_workspace

    SolidAngle(InputWorkspace=input_workspace,
               OutputWorkspace=solid_angle_ws,
               Method=detector_type)
    Divide(LHSWorkspace=input_workspace,
           RHSWorkspace=solid_angle_ws,
           OutputWorkspace=output_workspace)
    DeleteWorkspace(solid_angle_ws)
    ReplaceSpecialValues(InputWorkspace=output_workspace,
                         OutputWorkspace=output_workspace, NaNValue=0.,
                         InfinityValue=0.)

    return mtd[output_workspace]
