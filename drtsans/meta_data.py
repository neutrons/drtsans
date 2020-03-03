# Method in this module is to set meta data to SANS Mantid Workspaces
from drtsans.geometry import search_sample_detector_distance_meta_name


__all__ = ['set_meta_data', 'set_up_sample_detector_distance']


def set_up_sample_detector_distance(workspace, sample_detector_distance, distance_unit,
                                    non_exist_default_name):
    """Set detector to sample distance to meta data

    Parameters
    ----------
    workspace: str, ~mantid.api.MatrixWorkspace
        Mantid workspace instance or workspace name
    sample_detector_distance:  float, None
        sample to detector distance in meter
    distance_unit: str
        unit for source sample distance
    non_exist_default_name: str
        meta data name for sample detector distance if None of the default name exist

    Returns
    -------
    ~list
        item = (str, float, str) as log name, log value, log unit

    """
    # Ignore if source to sample distance is None
    if sample_detector_distance is None:
        return

    # Get the possible meta data names for this distance
    log_names = zip(*search_sample_detector_distance_meta_name(workspace, None))[0]

    # Set meta data name from default if the default one cannot be found
    if len(log_names) == 0:
        log_names = [non_exist_default_name]

    # Create the list for log, value, unit
    meta_overwrite_list = list()
    for log_name in log_names:
        meta_overwrite_list.append((log_name, sample_detector_distance, distance_unit))

    return meta_overwrite_list


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
        sample aperture radius in unit mm
    sample_thickness: None, float
        sample thickness in unit ???
    source_aperture_size: float, None
        source aperture size radius in unit mm
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter

    Returns
    -------

    """
    raise RuntimeError('set meta data is virtual in module drtsans')
