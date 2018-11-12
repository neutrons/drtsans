from mantid.simpleapi import ApplyTransmissionCorrection
from mantid.kernel import logger


def apply_transmission(input_ws, output_ws, trans_value=None, trans_error=None,
                       trans_ws=None, theta_dependent=true):
    '''
    Apply a transmission correction to 2D SANS data.

    Either use trans_value and trans_error or trans_ws

    input_ws - Workspace to apply the transmission correction to
    output_ws - Workspace to store the corrected data in
    trans_ws - Workspace containing the transmission values [optional]
    trans_value - Transmission value to apply to all wavelengths. If specified.
    trans_error - The error on the transmission value (default 0.0)
    theta_dependent - If true, a theta-dependent transmission correction will 
                      be applied.

    '''

    if trans_value is not None and trans_error is not None:
        ApplyTransmissionCorrection(
            InputWorkspace=input_ws,
            OutputWorkspace=output_ws,
            TransmissionValue=trans_value,
            TransmissionError=trans_error,
            ThetaDependent=theta_dependent,
        )
    elif trans_ws is not None:
        ApplyTransmissionCorrection(
            InputWorkspace=input_ws,
            OutputWorkspace=output_ws,
            TransmissionWorkspace=trans_ws,
            ThetaDependent=theta_dependent
        )
    else:
        logger.error("Input not valid: Use trans_value + trans_value"
                     " or trans_ws.")
