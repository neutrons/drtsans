# Method in this module is to set meta data to SANS Mantid Workspaces
from mantid.simpleapi import AddSampleLogMultiple, AddSampleLog


__all__ = ['set_meta_data']


def set_meta_data(workspace, wave_length=None, wavelength_spread=None,
                  sample_offset=0.,
                  sample_aperture_diameter=None, sample_thickness=None,
                  source_aperture_diameter=None,
                  smearing_pixel_size_x=None, smearing_pixel_size_y=None):
    """Set meta data to SANS Mantid Workspace as run properties

    Parameters
    ----------
    workspace: str, ~mantid.api.MatrixWorkspace
        Mantid workspace instance or workspace name
    wave_length: float, None
        wave length in Angstrom
    wavelength_spread: float, None
        wave length spread in Angstrom
    sample_offset: float
        offset of sample from origin in unit mm
    sample_aperture_diameter: float, None
        sample aperture diameter in mm
    sample_thickness: None, float
        sample thickness in unit cm
    source_aperture_diameter: float, None
        source aperture size radius in unit mm
    smearing_pixel_size_x: float, None
        pixel size in x direction in unit as meter, only for Q-resolution calculation
    smearing_pixel_size_y: float, None
        pixel size in Y direction in unit as meter, only for Q-resolution calculation

    Returns
    -------

    """
    # Init list for sample log name, value and unit
    meta_data_list = list()

    # Wave length and wave length spread shall be set Number Series
    # Wave length
    if wave_length is not None:
        # meta_data_list.append(('wavelength', np.array([wave_length, wave_length]), 'A'))
        AddSampleLog(workspace, LogName='wavelength', LogText='{}'.format(wave_length), LogType='Number Series',
                     LogUnit='A')

    # Wave length spread
    if wavelength_spread is not None:
        # meta_data_list.append(('wavelength_spread', np.array([wavelength_spread, wavelength_spread]), 'A'))
        AddSampleLog(workspace, LogName='wavelength_spread', LogText='{}'.format(wavelength_spread),
                     LogType='Number Series')

    # Add the sample log dictionary to add
    if sample_aperture_diameter is not None:
        meta_data_list.append(('sample_aperture_diameter', sample_aperture_diameter, 'mm'))

    # Source aperture radius
    if source_aperture_diameter is not None:
        meta_data_list.append(('source_aperture_diameter', source_aperture_diameter, 'mm'))

    # Sample offset
    meta_data_list.append(('sample_offset', sample_offset, 'mm'))

    # Sample thickness
    if sample_thickness is not None:
        meta_data_list.append(('sample_thickness', sample_thickness, 'cm'))

    # Pixel size
    if smearing_pixel_size_x is not None and smearing_pixel_size_y is not None:
        meta_data_list.append(('smearingPixelSizeX', smearing_pixel_size_x, 'm'))
        meta_data_list.append(('smearingPixelSizeY', smearing_pixel_size_y, 'm'))
    elif smearing_pixel_size_x is None and smearing_pixel_size_y is None:
        pass
    else:
        raise RuntimeError('Pixel size X ({}) and Y ({}) must be set together'
                           ''.format(smearing_pixel_size_x, smearing_pixel_size_y))

    # Add log value
    if len(meta_data_list) > 0:
        # only work on non-empty meta data list
        log_names, log_values, log_units = zip(*meta_data_list)

        # add meta data (as sample logs) to workspace
        AddSampleLogMultiple(Workspace=workspace, LogNames=log_names,
                             LogValues=log_values,
                             LogUnits=log_units)
