from mantid import mtd, logger

from drtsans.samplelogs import SampleLogs

__all__ = ['beam_radius',]


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
        Units of the output beam radius. Either 'mm' or 'm'.
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
    """
    from drtsans.tof.eqsans.geometry import beam_radius as eqsans_beam_radius
    ws = mtd[str(input_workspace)]
    if ws.getInstrument().getName() == 'EQ-SANS':
        return eqsans_beam_radius(ws, unit='mm')

    # Apertures, assumed to be in mili-meters
    sample_logs = SampleLogs(ws)
    radius_sample_aperture = sample_logs[sample_aperture_diameter_log].value / 2.
    radius_source_aperture = sample_logs[source_aperture_diameter_log].value / 2.

    # Distances
    ssd = sample_logs[ssd_log].value
    sdd = sample_logs[sdd_log].value

    # Calculate beam radius
    radius = radius_sample_aperture + sdd * (radius_sample_aperture + radius_source_aperture) / ssd

    logger.notice("Radius calculated from the input workspace = {:.2} mm".format(radius))
    return radius if unit == 'mm' else 1e-3 * radius