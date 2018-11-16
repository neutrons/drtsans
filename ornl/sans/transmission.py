from mantid.simpleapi import (
    ApplyTransmissionCorrection, FindDetectorsInShape, GroupDetectors, Divide,
    DeleteWorkspaces)
from mantid.kernel import logger


def apply_transmission(input_ws, output_ws, trans_value=None, trans_error=None,
                       trans_ws=None, theta_dependent=True):
    """
    Correct intensities with a transmission coefficient.

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace to apply the transmission correction to
    output_ws: str
        Name of the workspace to store the corrected data in
    trans_value: float
        Zero-angle transmission value to apply to all wavelengths.
    trans_error: float
        The error on the zero-angle transmission value (default 0.0)
    trans_ws: MatrixWorkspace
        Workspace containing the fitted or raw zero-angle transmission
        values [optional]
    theta_dependent: bool
        if True, a theta-dependent transmission correction will be applied.
    """

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


def _calculate_radius_from_input_ws(input_ws):
    """
    Calculate the radius according to:
    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Input workspace

    Returns
    -------
    float
        Radius
    """
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
    """
    Find the detectors within a circle of the beam center

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace containing the detector already beam-centered
    radius: float
        Radius of the circle encompassing the detectors of interest. Units
        in mili meters

    Returns
    -------
    numpy.ndarray
        Detector ID's
    """

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


def zero_angle_transmission(input_sample_ws, input_reference_ws,
                            output_ws, radius=None, delete_temp_wss=True):
    """
    Calculate the transmission coefficients at zero scattering angle

    Parameters
    ----------
    input_sample_ws: MatrixWorkspace
        Sample workspace (possibly obtained with an attenuated beam)
    input_reference_ws: MatrixWorkspace
    output_ws: str
        Name of the output workspace containing the transmission values.
    radius: float
        Radius around the bean center for pixel integration. If None,
        the beam radius is used, calculated using the sample
        and source apertures
    delete_temp_wss: bool
        Delete the grouping detector workspaces

    Returns
    -------
    MatrixWorkspace
        Workspace containing the raw transmission values (its name is
        given by `output_ws`)
    """
    if radius is None:
        radius = _calculate_radius_from_input_ws(input_reference_ws)

    detector_ids = _get_detector_ids_from_radius(input_reference_ws, radius)

    # by default it sums all the grouped detectors
    input_sample_ws_grouped = GroupDetectors(InputWorkspace=input_sample_ws,
                                             DetectorList=detector_ids)

    input_reference_ws_grouped = GroupDetectors(
        InputWorkspace=input_reference_ws,
        DetectorList=detector_ids)

    # Raw transmission values
    ws = Divide(LHSWorkspace=input_sample_ws_grouped,
                RHSWorkspace=input_reference_ws_grouped,
                OutputWorkspace=output_ws)

    if delete_temp_wss:
        DeleteWorkspaces(
            WorkspaceList=[input_sample_ws_grouped,
                           input_reference_ws_grouped])
    return ws
