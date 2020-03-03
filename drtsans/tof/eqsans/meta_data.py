# Method in this module is to set meta data to EQSANS Mantid Workspaces
from mantid.simpleapi import AddSampleLogMultiple
from drtsans.meta_data import set_up_sample_detector_distance, set_up_source_sample_distance

__all__ = ['set_meta_data']


def set_meta_data(workspace, wave_length=None, wavelength_spread=None,
                  sample_to_detector_distance=None, source_to_sample_distance=None,
                  sample_aperture_diameter=None, sample_thickness=None,
                  source_aperture_diameter=None,
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
        source to sample distance in mm
    sample_aperture_diameter: float, None
        sample aperture diameter in unit mm
    sample_thickness: None, float
        sample thickness in unit cm
    source_aperture_diameter: float, None
        source aperture diameter in unit meter
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter


    Returns
    -------

    """
    # Exception
    if wave_length is not None or wavelength_spread is not None:
        raise RuntimeError('Wave length and wave length spread are not allowed to set to EQ-SANS')

    # Log value dictionary: 3-tuple (log name, log value, unit)
    meta_data_list = list()

    # Add the sample log dictionary to add
    if sample_aperture_diameter is not None:
        meta_data_list.append(('sample_aperture_diameter', sample_aperture_diameter, 'mm'))

    # Source aperture radius
    if source_aperture_diameter is not None:
        meta_data_list.append(('source_aperture_diameter', source_aperture_diameter, 'mm'))

    # Source sample distance
    if source_to_sample_distance is not None:
        # meta_data_list.append(('source_sample_distance', source_to_sample_distance, 'm'))
        sub_list = set_up_source_sample_distance(workspace,  source_to_sample_distance, 'm',
                                                 'source_sample_distance')
        meta_data_list.extend(sub_list)

    # Sample to detector distance
    if sample_to_detector_distance is not None:
        sub_list = set_up_sample_detector_distance(workspace, source_to_sample_distance, 'm',
                                                   'detectorZ')
        meta_data_list.extend(sub_list)

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
