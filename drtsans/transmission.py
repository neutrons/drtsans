import numpy as np

from mantid import mtd
from mantid.kernel import logger
from mantid.simpleapi import (ApplyTransmissionCorrection, Divide, FindDetectorsInShape, GroupDetectors,
                              ReplaceSpecialValues, RebinToWorkspace)
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import masked_detectors


# To-do. This should be substituted with a function similar to
# eqsans.transmission.beam_radius so that
# eqsans.transmission.beam_radius is not required anymore
def beam_radius(input_workspace, unit='mm',
                sample_aperture_diameter_log='sample-aperture-diameter',
                source_aperture_diameter_log='source-aperture-diameter',
                sdd_log='sample-detector-distance',
                ssd_log='source-sample-distance'):
    """
    Calculate the radius in mm according to:
    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD

    Parameters
    ----------
    input_workspace: MatrixWorkspace, str
        Input workspace
    unit: str
        Either 'mm' or 'm'.
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
    from drtsans.tof.eqsans import beam_radius as eqsans_beam_radius
    ws = mtd[str(input_workspace)]
    if ws.getInstrument().getName() == 'EQ-SANS':
        return eqsans_beam_radius(ws, unit='mm')
    try:
        # Apertures, assumed to be in mili-meters
        sample_logs = SampleLogs(ws)
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
    return radius if unit == 'mm' else 1e-3 * radius


def detector_ids(input_workspace, radius, unit='mm'):
    """
    Find the detectors ID's within a certain radius from the beam center

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Workspace containing the detector already beam-centered
    radius: float
        Radius of the circle encompassing the detectors of interest.
    unit: str
        Either 'mm' or 'm', unit of the `radius` option.

    Returns
    -------
    numpy.ndarray
        List of detector ID's
    """
    r = radius * 1e-3 if unit == 'mm' else radius

    cylinder = r"""
    <infinite-cylinder id="shape">
        <centre x="0.0" y="0.0" z="0.0" />
        <axis x="0.0" y="0.0" z="1" />
        <radius val="{}" />
    </infinite-cylinder>
    <algebra val="shape" />
    """.format(r)
    det_ids = FindDetectorsInShape(Workspace=input_workspace,
                                   ShapeXML=cylinder)
    return det_ids


def calculate_transmission(input_sample, input_reference,
                           radius=None, radius_unit='mm',
                           output_workspace=None):
    """
    Calculate the raw transmission coefficients at zero scattering angle
    from already prepared sample and reference data.

    For EQ-SANS, one additional step fitting the returned raw values is
    necessary. Use `eqsans.calculate_transmission` instead.

    Parameters
    ----------
    input_sample: MatrixWorkspace
        Prepared sample workspace (possibly obtained with an attenuated beam)
    input_reference: MatrixWorkspace
        Prepared direct beam workspace (possibly obtained with an attenuated
         beam)
    radius: float
        Radius around the bean center for pixel integration, in millimeters.
        If None, radius will be obtained or calculated using `input_reference`.
    radius_unit: str
        Either 'mm' or 'm', and only used in conjunction with option `radius`.
    output_workspace: str
        Name of the output workspace containing the raw transmission values.
        If None, a hidden random name will be provided.

    Returns
    -------
    MatrixWorkspace
        Workspace containing the raw transmission values
    """
    if output_workspace is None:
        output_workspace = uwd()
    if radius is None:
        r = beam_radius(input_reference, unit='mm')
    else:
        r = radius if radius_unit == 'mm' else 1.e3 * radius  # to mm
    det_ids = detector_ids(input_reference, r, unit='mm')

    # Avoid pitfall of masking many pixels around the beam center
    msg = 'More than half of the detectors ' +\
          'within a radius of {:.2f} mm '.format(r) +\
          'from the beam center are masked in the input {0}'
    for k, v in dict(sample=input_sample, reference=input_reference).items():
        if len(masked_detectors(v, det_ids)) > len(det_ids) / 2:
            raise RuntimeError(msg.format(k))

    # by default it sums all the grouped detectors
    gis = GroupDetectors(InputWorkspace=input_sample, DetectorList=det_ids,
                         OutputWorkspace=uwd())
    gir = GroupDetectors(InputWorkspace=input_reference, DetectorList=det_ids,
                         OutputWorkspace=uwd())
    gir = RebinToWorkspace(WorkspaceToRebin=gir, WorkspaceToMatch=gis,
                           OutputWorkspace=gir.name())
    # calculate zero angle transmission coefficient(s)
    zat = Divide(LHSWorkspace=gis, RHSWorkspace=gir,
                 OutputWorkspace=output_workspace)
    av_t, av_e = np.mean(zat.dataY(0)), np.linalg.norm(zat.dataE(0))
    message = 'Average zero angle transmission = {0} +/- {1}'
    logger.notice(message.format(av_t, av_e))
    gis.delete()
    gir.delete()
    return zat


def apply_transmission_correction(input_workspace, trans_workspace=None, trans_value=None, trans_error=0.0,
                                  theta_dependent=True, output_workspace=None):
    r"""
    Correct the intensities with transmission coefficient(s).

    Parameters
    ----------
    input_workspace: MatrixWorkspace, str
        Input workspace to correct its intensities
    trans_workspace: MatrixWorkspace, str
        Workspace containing the transmission coefficient(s). The result of applying `calculate_transmission`
        to the input workspace. If None, `trans_value` will be used.
    trans_value: float
        A single transmission coefficient to correct the intensities. If None, `trans_workspace` will be used.
    trans_error: float
        Error associated to `trans_value`.
    output_workspace: str
        Name of the workspace containing the corrected intensities. If None, the `input_workspace`
        will be overwritten.

    Returns
    -------
    MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    kwargs = dict(InputWorkspace=input_workspace, ThetaDependent=theta_dependent, OutputWorkspace=output_workspace)
    if trans_workspace is not None:
        # EQ-SANS transmissions in skip-frame mode have transmission values of zero in the wavelength gap.
        # Need to be replaced with one to avoid division of intensities by zero.
        clean_trans_workspace = ReplaceSpecialValues(InputWorkspace=trans_workspace, SmallNumberThreshold=1.0e-6,
                                                     SmallNumberValue=1.0, OutputWorkspace=uwd())
        kwargs['TransmissionWorkspace'] = clean_trans_workspace
    elif trans_value is not None:
        kwargs.update(dict(TransmissionValue=trans_value, TransmissionError=trans_error))
    else:
        raise RuntimeError('Provide either trans_workspace or trans_value')
    ApplyTransmissionCorrection(**kwargs)
    if trans_workspace is not None:
        clean_trans_workspace.delete()
    return mtd[output_workspace]
