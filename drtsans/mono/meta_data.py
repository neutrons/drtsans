# Method in this module is to set meta data to SANS Mantid Workspaces
from mantid.simpleapi import AddSampleLogMultiple
from mantid.simpleapi import mtd


__all__ = ['set_meta_data']


def set_meta_data(workspace, wave_length=None, wavelength_spread=None,
                  sample_offset=0.,
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
    sample_offset: float
        offset of sample from origin in unit mm
    sample_aperture_diameter: float, None
        sample aperture diameter in mm
    sample_thickness: None, float
        sample thickness in unit cm
    source_aperture_diameter: float, None
        source aperture size radius in unit mm
    pixel_size_x: float, None
        pixel size in x direction in unit as meter
    pixel_size_y: float, None
        pixel size in Y direction in unit as meter

    Returns
    -------

    """
    import numpy as np

    # Init list for sample log name, value and unit
    meta_data_list = list()

    # Wave length
    if wave_length is not None:
        meta_data_list.append(('wavelength', np.array([wave_length, wave_length]), 'A'))

    # Wave length spread
    if wavelength_spread is not None:
        meta_data_list.append(('wavelength_spread', np.array([wavelength_spread, wavelength_spread]), 'A'))

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
    if pixel_size_x is not None and pixel_size_y is not None:
        meta_data_list.append(('pixel_size_x', pixel_size_x, 'm'))
        meta_data_list.append(('pixel_size_y', pixel_size_y, 'm'))
    elif pixel_size_x is None and pixel_size_y is None:
        pass
    else:
        raise RuntimeError('Pixel size X ({}) and Y ({}) must be set together'
                           ''.format(pixel_size_x, pixel_size_y))

    # Add log value
    if len(meta_data_list) > 0:
        # only work on non-empty meta data list
        log_names, log_values, log_units = zip(*meta_data_list)

        #  wavelength = runObj['wavelength'].getStatistics().mean
        #  wavelength_spread = runObj['wavelength_spread'].getStatistics().mean

        print('DEBUG wzz: Overwriting {} with value {}'
              ''.format(log_names, log_values))

        # add meta data (as sample logs) to workspace
        AddSampleLogMultiple(Workspace=workspace, LogNames=log_names,
                             LogValues=log_values,
                             LogUnits=log_units)

        # Check
        output_ws = mtd[str(workspace)]
        ws_run = output_ws.getRun()
        print(ws_run['wavelength'])
