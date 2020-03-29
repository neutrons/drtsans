# Method in this module is to set meta data to SANS Mantid Workspaces
from mantid.simpleapi import AddSampleLogMultiple, AddSampleLog
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import sample_detector_distance

__all__ = ['set_meta_data', 'get_sample_detector_offset']


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

        # add meta data (as sample logs) to workspace
        AddSampleLogMultiple(Workspace=workspace, LogNames=log_names,
                             LogValues=log_values,
                             LogUnits=log_units)


def get_sample_detector_offset(workspace, sample_si_meta_name, zero_sample_offset_sample_si_distance,
                               overwrite_sample_si_distance=None, overwrite_sample_detector_distance=None):
    """Get sample offset and detector offset from meta data

    This method is based on the assumption and fact that
    "sample position is set to nominal position (0, 0, 0) regardless of sample log SampleToSi"
    "detector1 is set to [0, 0, sample_detector_distance]
    It will be re-implemented

    Parameters
    ----------
    workspace: str, ~mantid.api.MatrixWorkspace
        Mantid workspace instance or workspace name
    sample_si_meta_name : str
        Sample to Si (window) meta data name
    zero_sample_offset_sample_si_distance: float
        default sample to Si window distance, i.e., distance without sample offset. unit = meter
    overwrite_sample_si_distance: float or None
        sample to Si window distance to overwrite.  Unit = mm (consistent with the unit of original meta data)
    overwrite_sample_detector_distance : float or None
        sample detector distance to overwrite. Unit = m (consistent with the unit of original meta data)

    Returns
    -------
    ~tuple
        sample offset (float) in unit meter and detector offset (float) in unit meter

    """
    # Calculate the sample offset and detector offset without overwriting value
    # This is caused by incorrect IDF which does not use SampleToSi.
    sample_logs = SampleLogs(workspace)
    # read sample log for SampleToSi and convert to meter from mm
    sample_to_si = sample_logs.find_log_with_units(sample_si_meta_name, 'mm') * 1E-3
    print('[DEBUG] Meta-Data Sample to Si = {} meter'.format(sample_to_si))
    print('[DEBUG] Hardcoded Sample to nominal distance = {} meter'.format(zero_sample_offset_sample_si_distance))

    # Offsets: shift both sample and detector to conserve sample-detector distance
    # Move instrument_component sample (relative) to [0, 0, 0.071 - SampleToSi/1000]
    sample_offset = zero_sample_offset_sample_si_distance - sample_to_si
    # Move instrument_component detector1 relative [0, 0, 0.071 - SampleToSi/1000]
    detector_offset = sample_offset

    # Get sample detector distance by calculation from instrument geometry directly
    sample_det_distance = sample_detector_distance(workspace, unit='m')
    print('[DEBUG] Calculated sample detector distance = {}'.format(sample_det_distance))

    # With overwriting distance(s)
    if overwrite_sample_si_distance is not None or overwrite_sample_detector_distance is not None:
        # 2 cases to handle.  The order must be conserved
        if overwrite_sample_si_distance is not None:
            # Sample-Si distance is overwritten. NeXus-recorded sample-detector-distance is thus inaccurate.
            # # convert unit of (overwrite)-sample-Si-distance to meter
            # overwrite_sample_si_distance *= 1E-3

            # Shift the sample position only without moving detector
            overwrite_offset = sample_to_si - overwrite_sample_si_distance
            print('[DEBUG]: SampleToSi = {}, SampleToSiOverwrite = {}, Original SampleOffset = {}'
                  ''.format(sample_to_si, overwrite_sample_si_distance, sample_offset))
            sample_offset += overwrite_offset
            sample_det_distance -= overwrite_offset
            sample_to_si = overwrite_offset

        if overwrite_sample_detector_distance is not None:
            # Sample-detector distance is overwritten, i.e., fix the sample position and move detector to
            # make the sample-detector-distance to the overwriting value
            # Move instrument_component detector1 relatively
            # [0, 0, sample_detector_distance_overwrite - sample_detector_distance_nexus]
            overwrite_offset = overwrite_sample_detector_distance - sample_det_distance
            detector_offset += overwrite_offset
            sample_det_distance += overwrite_offset
    # END-IF

    print('[INFO] Sample offset = {}, Detector offset = {}, Sample-detector-distance = {}, '
          'Sample-Si-window-distance = {}'
          ''.format(sample_offset, detector_offset, sample_det_distance, sample_to_si))

    return sample_offset, detector_offset
