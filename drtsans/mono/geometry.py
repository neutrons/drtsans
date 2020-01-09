from mantid import logger
from drtsans.geometry import (sample_aperture_diameter, sample_detector_distance, source_aperture_diameter,
                              source_sample_distance)

__all__ = ['beam_radius']


def beam_radius(input_workspace, unit='mm'):
    """
    Calculate the beam radius impinging on the detector

    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD, where
    R_sampleAp: radius of the sample aperture,
    SDD: distance between the sample and the detector,
    R_sourceAp: radius of the source aperture,
    SSD: distance between the source and the sample.

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace, str
        Input workspace
    unit: str
        Units of the output beam radius. Either 'mm' or 'm'.

    Returns
    -------
    float
    """
    r_sa = sample_aperture_diameter(input_workspace, unit='m') / 2.0  # radius
    r_so = source_aperture_diameter(input_workspace, unit='m') / 2.0  # radius
    l1 = source_sample_distance(input_workspace, unit='m')
    l2 = sample_detector_distance(input_workspace, unit='m')

    radius = r_sa + (r_sa + r_so) * (l2 / l1)

    logger.notice("Radius calculated from the input workspace = {:.2} mm".format(radius * 1e3))
    return radius if unit == 'm' else 1e3 * radius
