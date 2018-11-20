from mantid.simpleapi import (
    ApplyTransmissionCorrection, FindDetectorsInShape, GroupDetectors, Divide,
    DeleteWorkspaces)
from mantid.kernel import logger
from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogs


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

    Returns
    -------
    namedtuple or None
        Return value of Mantid's algorithm ApplyTransmissionCorrection. None
        if input is invalid
    """

    if trans_value is not None and trans_error is not None:
        return ApplyTransmissionCorrection(InputWorkspace=input_ws,
                                           OutputWorkspace=output_ws,
                                           TransmissionValue=trans_value,
                                           TransmissionError=trans_error,
                                           ThetaDependent=theta_dependent)
    elif trans_ws is not None:
        return ApplyTransmissionCorrection(InputWorkspace=input_ws,
                                           OutputWorkspace=output_ws,
                                           TransmissionWorkspace=trans_ws,
                                           ThetaDependent=theta_dependent)
    else:
        logger.error("Input not valid: Use trans_value + trans_value"
                     " or trans_ws.")
        return None


def calculate_radius_from_input_ws(input_ws,
                                   sample_ad='sample_aperture-diameter',
                                   source_ad='source_aperture-diameter',
                                   sdd_log='sample-detector-distance',
                                   sasd_log='source-sample-distance'):
    """
    Calculate the radius according to:
    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Input workspace
    sample_ad: str
        Log entry for the sample-aperture diameter
    source_ad: str
        Log entry for the source-aperture diameter
    sdd_log: str
        Log entry for the sample to detector distance
    sasd_log: str
        Log entry for the source-aperture to sample distance

    Returns
    -------
    float
        Radius, in mili-meters
    """
    try:
        # Apertures
        sl = SampleLogs(input_ws)
        radius_sample_aperture = sl[sample_ad].value / 2.
        radius_source_aperture = sl[source_ad].value / 2.
        # Distances
        sasd = sl[sasd_log].value
        sdd = sl[sdd_log].value
        # Calculate beam radius
        radius = radius_sample_aperture + sdd * \
            (radius_sample_aperture + radius_source_aperture) / sasd
    except ValueError as error:
        logger.error(
            "Some of the properties are likely to not exist in the WS")
        logger.error(error)
        raise
    logger.notice("Radius calculated from the WS = {}".format(radius))
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


@namedtuplefy
def zero_angle_transmission(input_sample_ws, input_reference_ws, radius,
                            output_ws, delete_temp_wss=True):
    """
    Calculate the raw transmission coefficients at zero scattering angle

    Parameters
    ----------
    input_sample_ws: MatrixWorkspace
        Sample workspace (possibly obtained with an attenuated beam)
    input_reference_ws: MatrixWorkspace
        Direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        Radius around the bean center for pixel integration, in mili-meters
    output_ws: str
        Name of the output workspace containing the transmission values.
    delete_temp_wss: bool
        Delete the grouping detector workspaces

    Returns
    -------
    namedtuple
        Fields of the namedtuple
        - transmission: MatrixWorkspace, Workspace containing the raw transmission
            values (its name is given by `output_ws`)
        - radius: float, integration radius, in pixel units
        - detids: list, list of detector ID's used for integration
        - sample: MatrixWorkspace, sample workspace after integration.
            `None` if `delete_temp_wss` is True
        - reference: MatrixWorkspace, reference workspace after integration.
            `None` if `delete_temp_wss` is True
    """

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
        input_sample_ws_grouped = None
        input_reference_ws_grouped = None

    return dict(transmission=ws,
                detids=list(detector_ids),
                sample=input_sample_ws_grouped,
                reference=input_reference_ws_grouped)
