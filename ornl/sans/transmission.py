from mantid.simpleapi import (
    ApplyTransmissionCorrection, FindDetectorsInShape, GroupDetectors, Divide,
    DeleteWorkspaces)
from mantid.kernel import logger


def apply_transmission(input_ws, output_ws, trans_value=None, trans_error=None,
                       trans_ws=None, theta_dependent=True):
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
                    


#
#
#

def _calculate_radius_from_input_ws(input_ws):
    '''
    Calculate the radius according to:
    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD
    '''

    r = input_ws.getRun()

    try:
        radius_sample_aperture = r.getProperty(
            "sample-aperture-diameter").value
        radius_source_aperture = r.getProperty(
            "source-aperture-diameter").value
        ssd = r.getProperty("source-sample-distance").value
        sdd = r.getProperty("sample-detector-distance").value

        radius = radius_sample_aperture + sdd * \
            (radius_sample_aperture + radius_source_aperture) / ssd
    except ValueError as error:
        logger.error(
            "Some of the properties are likely to not exist in the WS")
        logger.error(error)
        raise
    logger.notice("Radis calculated from the WS = {}".format(radius))
    return radius


def _get_detector_ids_from_radius(input_ws, radius):
    '''
    Radius is in mm
    '''
    radius_in_meters = radius * 1e-3
    cylinder = """
    <infinite-cylinder id=\"shape\">
        <centre x=\"0.0\" y=\"0.0\" z=\"0.0\" />
        <axis x=\"0.0\" y=\"0.0\" z=\"1\" />
        <radius val=\"{}\" />
    </infinite-cylinder>
    <algebra val=\"shape\" />
    """.format(radius_in_meters)

    detector_ids = FindDetectorsInShape(input_ws, cylinder)
    return detector_ids


def calculate_transmission(input_sample_ws, input_reference_ws, output_ws,
                           radius=None, delete_temp_wss=True):
    '''
    If the radius is none calculates it according to
    _calculate_radius_from_input_ws.

    Creates a transmission Workspace: output_ws
    '''

    if radius is None:
        radius = _calculate_radius_from_input_ws(input_reference_ws)

    detector_ids = _get_detector_ids_from_radius(input_reference_ws, radius)

    # by default it sums all the grouped detectors
    input_sample_ws_grouped = GroupDetectors(InputWorkspace=input_sample_ws,
                                             DetectorList=detector_ids)

    input_reference_ws_grouped = GroupDetectors(
        InputWorkspace=input_reference_ws,
        DetectorList=detector_ids)

    output_ws = Divide(LHSWorkspace=input_sample_ws_grouped,
                       RHSWorkspace=input_reference_ws_grouped)

    if delete_temp_wss:
        DeleteWorkspaces(
            WorkspaceList=[input_sample_ws_grouped,
                           input_reference_ws_grouped])
    return output_ws
