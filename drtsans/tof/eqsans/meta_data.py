# Method in this module is to set meta data to EQSANS Mantid Workspaces
from mantid.api import AnalysisDataService as mtd
from mantid.simpleapi import AddSampleLogMultiple

__all__ = ['set_meta_data']


def set_meta_data(workspace, wave_length=None, wavelength_spread=None,
                  sample_to_detector_distance=None, source_to_sample_distance=None,
                  sample_aperture_size=None, sample_thickness=None,
                  source_aperture_size=None,
                  pixel_size_x=None, pixel_size_y=None):
    """Set meta data to SANS Mantid Workspace as run properties

    Parameters
    ----------
    workspace: str, ~mantid.api.MatrixWorkspace
        Mantid workspace instance or workspace name
    wave_length: float, None
        wave length in Angstrom
    wavelength_spread: float, None
        wave length spread in Angstrom
    sample_to_detector_distance: float, None
        sample to detector distance in meter
    source_to_sample_distance: float, None
        source to sample distance in meter
    sample_aperture_size: float, None
        sample aperture size (radius or diameter????)
    sample_thickness: None, float
        sample thickness in unit ???
    source_aperture_size: float, None
        source aperture size (radius ??? diameter???) in unit ????
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter


    Returns
    -------

    """
    # Get workspace
    if isinstance(workspace, str):
        workspace = mtd[workspace]

    # Exception
    if wave_length is None or wavelength_spread is not None:
        raise RuntimeError('Wave length and wave length spread are not allowed to set to EQ-SANS')

    # Log value dictionary: 3-tuple (log name, log value, unit)
    meta_data_list = list()

    # Add the sample log dictionary to add
    if sample_aperture_size is not None:
        meta_data_list.append(('sample_aperture_radius', sample_aperture_size, 'mm'))

    # Source aperture radius
    if source_aperture_size is not None:
        meta_data_list.append(('source_aperture_radius', source_aperture_size, 'mm'))

    # Source sample distance
    if source_to_sample_distance is not None:
        meta_data_list.append(('source_sample-distance', source_to_sample_distance, 'm'))

    # Sample to detector distance
    if sample_to_detector_distance is not None:
        meta_data_list.append(('detectorZ', sample_to_detector_distance, 'm'))

    # Add log value
    log_names, log_values, log_units = zip(*meta_data_list)

    AddSampleLogMultiple(Workspace=workspace, LogNames=log_names,
                         LogValues=log_values,
                         LogUnits=log_units)
