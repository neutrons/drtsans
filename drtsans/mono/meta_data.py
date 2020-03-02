# Method in this module is to set meta data to SANS Mantid Workspaces
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
        sample thickness in unit cm
    source_aperture_size: float, None
        source aperture size (radius ??? diameter???) in unit ????
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter

    Returns
    -------

    """
    # Init list for sample log name, value and unit
    meta_data_list = list()

    # Wave length
    if wave_length is not None:
        meta_data_list.append(('wavelength', wave_length, 'A'))

    # Wave length spead
    if wavelength_spread is not None:
        meta_data_list.append(('wavelength_spread', wavelength_spread, 'A'))

    # Add the sample log dictionary to add
    if sample_aperture_size is not None:
        meta_data_list.append(('sample_aperture_radius', sample_aperture_size, 'mm'))

    # Source aperture radius
    if source_aperture_size is not None:
        meta_data_list.append(('source_aperture_radius', source_aperture_size, 'mm'))

    # Source sample distance
    if source_to_sample_distance is not None:
        meta_data_list.append(('source_sample_distance', source_to_sample_distance, 'm'))

    # Sample to detector distance
    if sample_to_detector_distance is not None:
        meta_data_list.append(('sample_detector_distance', sample_to_detector_distance, 'm'))

    # Sample thickness
    if sample_thickness is not None:
        meta_data_list.append(('sample_thickness', sample_thickness, 'cm'))

    # Pixel size
    if pixel_size_x is not None and pixel_size_y is not None:
        meta_data_list.append(('pixel_size_x', pixel_size_x, 'm'))
        meta_data_list.append(('pixel_size_y', pixel_size_y, 'm'))
    elif pixel_size_x is None and pixel_size_y is None:
        pass
    else:
        raise RuntimeError('Pixel size X ({}) and Y ({}) must be set together'
                           ''.format(pixel_size_x, pixel_size_y))

    # Add log value
    log_names, log_values, log_units = zip(*meta_data_list)

    AddSampleLogMultiple(Workspace=workspace, LogNames=log_names,
                         LogValues=log_values,
                         LogUnits=log_units)
