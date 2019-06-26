from __future__ import absolute_import, division, print_function

from mantid import mtd
from mantid.kernel import logger
from mantid.simpleapi import (ApplyTransmissionCorrection,
                              FindDetectorsInShape, GroupDetectors,
                              ReplaceSpecialValues)
from ornl.sans.samplelogs import SampleLogsReader


def _apply_transmission_mantid(input_ws, trans_value=None, trans_error=None,
                               trans_ws=None, theta_dependent=True):
    """
    Correct intensities with a transmission coefficient.
    This justs calls the Mantid algorithm

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

    # I don't know why we still have to specify the OutputWorkspace...
    trans_corrected_ws_name = "__trans_corrected_ws"
    if trans_value is not None and trans_error is not None:
        ApplyTransmissionCorrection(
            InputWorkspace=input_ws, OutputWorkspace=trans_corrected_ws_name,
            TransmissionValue=trans_value, TransmissionError=trans_error,
            ThetaDependent=theta_dependent)
        return mtd[trans_corrected_ws_name]
    elif trans_ws is not None:
        ApplyTransmissionCorrection(
            InputWorkspace=input_ws, OutputWorkspace=trans_corrected_ws_name,
            TransmissionWorkspace=trans_ws, ThetaDependent=theta_dependent)
        return mtd[trans_corrected_ws_name]
    else:
        logger.error("Input not valid: Use trans_value + trans_value"
                     " or trans_ws.")
        return None


def _calculate_radius_from_input_ws(
        input_ws, sample_aperture_diameter_log='sample-aperture-diameter',
        source_aperture_diameter_log='source-aperture-diameter',
        sdd_log='sample-detector-distance', ssd_log='source-sample-distance'):
    """
    Calculate the radius in mm according to:
    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Input workspace
    sample_aperture_diameter_log: str
        Log entry for the sample-aperture diameter
    source_aperture_diameter_log: str
        Log entry for the source-aperture diameter
    sdd_log: str
        Log entry for the sample to detector distance
    ssd_log: str
        Log entry for the source-aperture to sample distance

    Returns
    -------
    float
        Radius, in millimeters
    """
    try:
        # Apertures
        sample_logs = SampleLogsReader(input_ws)
        radius_sample_aperture = sample_logs[
            sample_aperture_diameter_log].value / 2.
        radius_source_aperture = sample_logs[
            source_aperture_diameter_log].value / 2.
        # Distances
        ssd = sample_logs[ssd_log].value
        sdd = sample_logs[sdd_log].value
        # Calculate beam radius
        radius = radius_sample_aperture + sdd * \
            (radius_sample_aperture + radius_source_aperture) / ssd
    except ValueError as error:
        logger.error(
            "Some of the properties are likely to not exist in the WS")
        logger.error(error)
        raise
    logger.notice("Radius calculated from the WS = {:.2} mm".format(radius))
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
        in milimeters

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


def _zero_angle_transmission(input_sample_ws, input_reference_ws, radius):
    """
    Calculate the raw transmission coefficients at zero scattering angle

    Parameters
    ----------
    input_sample_ws: MatrixWorkspace
        Sample workspace (possibly obtained with an attenuated beam)
    input_reference_ws: MatrixWorkspace
        Direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        Radius around the bean center for pixel integration, in millimeters

    Returns
    -------
    MatrixWorkspace
        Workspace containing the raw transmission values
    """

    detector_ids = _get_detector_ids_from_radius(input_reference_ws, radius)

    # by default it sums all the grouped detectors
    __input_sample_ws_grouped = GroupDetectors(
        InputWorkspace=input_sample_ws, DetectorList=detector_ids)

    __input_reference_ws_grouped = GroupDetectors(
        InputWorkspace=input_reference_ws, DetectorList=detector_ids)

    __zero_angle_transmission_ws = __input_sample_ws_grouped / \
        __input_reference_ws_grouped

    logger.notice("Zero Angle Tranmission = {} +/-{}".format(
        __zero_angle_transmission_ws.readY(0)[0],
        __zero_angle_transmission_ws.readE(0)[0]))

    return __zero_angle_transmission_ws


def calculate_transmission_ws(input_sample_ws, input_reference_ws, radius=None,
                              theta_dependent=False):
    '''Only calculates the transmission using a radius and the zero angle
    transmission algorithm

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace to apply the transmission correction to
    input_reference_ws: MatrixWorkspace
        Direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        radius of the beam in mm. Eg. 50 mm. If None calculates it from the
        reference workspace
    theta_dependent: bool
        if True, a theta-dependent transmission correction will be applied.

    Returns
    -------
    MatrixWorkspace
        Returns a Mantid tranmission workspace
    '''
    if not radius:
        radius = _calculate_radius_from_input_ws(input_reference_ws)

    # This returns a WS with the sensitivity value + error
    __calculated_transmission_ws = _zero_angle_transmission(
        input_sample_ws, input_reference_ws, radius)

    return __calculated_transmission_ws


def calculate_transmission_value(input_sample_ws, input_reference_ws,
                                 radius=None, theta_dependent=False):
    '''Only calculates the transmission using a radius and the zero angle
    transmission algorithm

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace to apply the transmission correction to
    input_reference_ws: MatrixWorkspace
        Direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        radius of the beam in mm. Eg. 50 mm. If None calculates it from the
        reference workspace
    theta_dependent: bool
        if True, a theta-dependent transmission correction will be applied.

    Returns
    -------
    Tuple
        Tuple with transmission value and error: (T, sigma(T))
    '''
    # This returns a WS with the sensitivity value + error
    __calculated_transmission_ws = calculate_transmission_ws(
        input_sample_ws, input_reference_ws, radius, theta_dependent)

    return (__calculated_transmission_ws.readY(0)[0],
            __calculated_transmission_ws.readE(0)[0])


def apply_transmission_correction_ws(input_sample_ws, transmission_ws,
                                     theta_dependent=False):
    '''
    This is the main method used to correct for transmission

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace to apply the transmission correction to
    transmission_ws: MatrixWorkspace
        The Workspace calculated by Mantid with the transmission value / Error.
        Note that for EQSANS those are spectra not values!
    theta_dependent: bool
        if True, a theta-dependent transmission correction will be applied.

    Returns
    -------
    MatrixWorkspace
        The data corrected for transmission.
    '''

    __input_sample_trans_corrected = _apply_transmission_mantid(
        input_sample_ws, trans_ws=transmission_ws,
        theta_dependent=theta_dependent)

    __input_sample_trans_corrected_no_nans = ReplaceSpecialValues(
        InputWorkspace=__input_sample_trans_corrected, NaNValue=0,
        InfinityValue=0)

    return __input_sample_trans_corrected_no_nans


def apply_transmission_correction_value(input_sample_ws, trans_value,
                                        trans_error, theta_dependent=False):
    '''This is the main method used to correct for transmission

    Parameters
    ----------
    input_ws: MatrixWorkspace
        Workspace to apply the transmission correction to
    trans_value: float
        Zero-angle transmission value to apply to all wavelengths.
    trans_error: float
        The error on the zero-angle transmission value
    theta_dependent: bool
        if True, a theta-dependent transmission correction will be applied.

    Returns
    -------
    MatrixWorkspace
        The data corrected for transmission.
    '''

    __input_sample_trans_corrected = _apply_transmission_mantid(
        input_sample_ws,  trans_value=trans_value, trans_error=trans_error,
        theta_dependent=theta_dependent)

    __input_sample_trans_corrected_no_nans = ReplaceSpecialValues(
        InputWorkspace=__input_sample_trans_corrected, NaNValue=0,
        InfinityValue=0)

    return __input_sample_trans_corrected_no_nans
