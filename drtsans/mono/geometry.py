from mantid import mtd, logger


__all__ = ['beam_radius']


def beam_radius(input_workspace, unit='mm',
                sample_aperture_diameter_log='sample_aperture_diameter',
                source_aperture_diameter_log='source_aperture_diameter',
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
        Log entry for the sample_aperture diameter
    source_aperture_diameter_log: str
        Log entry for the source_aperture diameter
    sdd_log: str
        Log entry for the sample to detector distance
    ssd_log: str
        Log entry for the source_aperture to sample distance

    Returns
    -------
    float
    """
    from drtsans.tof.eqsans.geometry import beam_radius as eqsans_beam_radius
    from drtsans.mono.momentum_transfer import retrieve_instrument_setup
    ws = mtd[str(input_workspace)]
    if ws.getInstrument().getName() == 'EQ-SANS':
        return eqsans_beam_radius(ws, unit='mm')

    inst = retrieve_instrument_setup(ws)
    radius = inst.sample_aperture_radius + inst.sample_det_center_distance * (inst.sample_aperture_radius
                                                                              + inst.source_aperture_radius) / inst.l1

    logger.notice("Radius calculated from the input workspace = {:.2} mm".format(radius * 1e3))
    return radius if unit == 'm' else 1e3 * radius
