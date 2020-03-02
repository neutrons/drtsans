# Method in this module is to set meta data to SANS Mantid Workspaces


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
