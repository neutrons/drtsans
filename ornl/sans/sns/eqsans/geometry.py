from mantid.simpleapi import MoveInstrumentComponent
from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (get_instrument, detector_name,
                                source_sample_distance)

__all__ = ['detector_z_log',
           'translate_sample_by_z', 'translate_detector_by_z',
           'sample_aperture_diameter', 'source_aperture_diameter']

detector_z_log = 'detectorZ'


def translate_sample_by_z(ws, z):
    r"""
    Shift the position of the sample by the desired amount

    Parameters
    ----------
    ws: ~mantid.api.MatrixWorkspace
        Input workspace containing instrument file
    z: float
        Translation to be applied
    """
    MoveInstrumentComponent(Workspace=str(ws), Z=z,
                            ComponentName='sample-position',
                            RelativePosition=True)


def translate_detector_z(ws, z=None, relative=True):
    r"""
    Adjust Z-coordinate of detector bank in instrument file.


    Parameters
    ----------
    ws: ~mantid.api.MatrixWorkspace
        Input workspace containing instrument file
    z: float
        Translation to be applied, in units of meters. If :py:obj:`None`, log_key
        stored in ``detector_z_log`` is used
    relative: bool
        If :py:obj:`True`, add to the current z-coordinate. If :py:obj:`False`, substitute
        the current z-coordinate with the new value.
    """
    if z is None:
        sl = SampleLogs(ws)
        z = 1e-3 * sl.single_value(detector_z_log)  # assumed in mili-meters

    MoveInstrumentComponent(Workspace=ws, Z=z, ComponentName=detector_name(ws),
                            RelativePosition=relative)


def translate_detector_by_z(ws, z, **kwargs):
    r"""
    Simplified call to :func:`.translate_detector_z`
    """
    return translate_detector_z(ws, z=z, relative=True, **kwargs)


def source_monitor_distance(source, unit='mm', log_key=None, search_logs=True):
    r"""
    Report the distance (always positive!) between source and monitor.

    If logs are not used or distance fails to be found in the logs, then
    calculate the distance using the instrument configuration file.

    Parameters
    ----------
    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number
    unit: str
        'mm' (millimeters), 'm' (meters)
    log_key: str
        Only search for the given string in the logs. Do not use default
        log keys
    search_logs: bool
        Report the value found in the logs.

    Returns
    -------
    float
        distance between source and sample, in selected units
    """
    m2units = dict(mm=1e3, m=1.0)
    mm2units = dict(mm=1.0, m=1e-3)
    sl = SampleLogs(source)

    # Search the logs for the distance
    if search_logs is True:
        lk = 'source-monitor-distance' if log_key is not None else log_key
        try:
            return sl.single_value(lk) * mm2units[unit]
        except Exception:
            pass

    # Calculate the distance using the instrument definition file
    instrument = get_instrument(source)
    monitor = instrument.getComponentByName('monitor1')
    smd = monitor.getDistance(instrument.getSource())

    # Insert in the logs if not present
    if 'source-monitor-distance' not in sl.keys():
        sl.insert('source-monitor-distance', smd * 1.e3, unit='mm')

    return smd * m2units[unit]


def sample_aperture_diameter(run, unit='m'):
    r"""
    Find the sample aperture diameter from the logs.

    Log keys searched are 'sample-aperture-diameter' (override beamslit4) and 'beamslit4'.

    Parameters
    ----------
    run: Mantid Run instance, :py:obj:`~mantid.api.MatrixWorkspace`, file name, run number
        Input from which to find the aperture
    unit: str
        return aperture in requested length unit, either 'm' or 'mm'

    Returns
    -------
    float
        Sample aperture diameter, in requested units
    """
    sl = SampleLogs(run)
    sad = None
    for log_key in ('sample-aperture-diameter', 'beamslit4'):
        if log_key in sl.keys():
            sad = sl.single_value(log_key)
            break
    if sad is None:
        pnames = [p.name for p in run.run().getProperties()]
        raise RuntimeError('Unable to retrieve sample aperture diameter as neither log "sample-aperture-diameter" '
                           'nor "beamslit4" is in the sample logs.  Available logs are {}'
                           ''.format(pnames))

    if 'sample-aperture-diameter' not in sl.keys():
        sl.insert('sample-aperture-diameter', sad, unit='mm')
    if unit == 'm':
        sad /= 1000.0
    return sad


@namedtuplefy
def source_aperture(other, unit='m'):
    r"""
    Find the source aperture diameter and position

    After the moderator (source) there are three consecutive discs
    (termed wheels), each with eight holes in them (eight slits).
    Appropriate log entries (VbeamSlit, VbeamSlit2, VbeamSlit3) indicate
    the slit index for each of the three wheels. Thus, the position
    of the source and the source aperture are not the same. The most
    restrictive slit will define the source aperture

    Log entries beamslit, beamslit2, and beamslit3 store the required
    rotation angle for each wheel in order to align the appropriate slit
    with the neutron beam. These angles are not used in reduction.

    Parameters
    ----------
    other: Run, MatrixWorkspace, file name, run number
    unit: str
        Length unit, either 'm' or 'mm'

    Returns
    -------
    namedtuple
        Fields of the name tuple
        - float: diameter, in requested units
        - float: distance to sample, in requested units
    """
    n_wheels = 3
    index_to_diameter = [[5.0, 10.0, 10.0, 15.0, 20.0, 20.0, 25.0, 40.0],
                         [0.0, 10.0, 10.0, 15.0, 15.0, 20.0, 20.0, 40.0],
                         [0.0, 10.0, 10.0, 15.0, 15.0, 20.0, 20.0, 40.0]]
    distance_to_source = [10080, 11156, 12150]  # in mili-meters
    sl = SampleLogs(other)
    slit_indexes = [int(sl[log_key].value.mean()) - 1 for log_key in
                    ['vBeamSlit', 'vBeamSlit2', 'vBeamSlit3']]
    diameter = 20.0  # default slit size
    asd = -1.0  # aperture to sample distance
    ssd = source_sample_distance(other, unit='mm')
    for wheel_index in range(n_wheels):
        slit_index = slit_indexes[wheel_index]
        y = ssd - distance_to_source[wheel_index]  # aperture to sample dist
        if 0 <= slit_index < 6:
            x = index_to_diameter[wheel_index][slit_index]
            if asd < 0 or x / y < diameter / asd:
                diameter = x
                asd = y
    if unit == 'm':
        diameter /= 1000.0
        asd /= 1000.0
    return dict(diameter=diameter, distance_to_sample=asd, unit=unit)


def source_aperture_diameter(run, unit='mm'):
    r"""
    Find the source aperture diameter

    Either report log vale or compute this quantity.

    After the moderator (source) there are three consecutive discs
    (termed wheels), each with eight holes in them (eight slits).
    Appropriate log entries (VbeamSlit, VbeamSlit2, VbeamSlit3) indicate
    the slit index for each of the three wheels. Thus, the position
    of the source and the source aperture are not the same. The most
    restrictive slit will define the source aperture

    Log entries beamslit, beamslit2, and beamslit3 store the required
    rotation angle for each wheel in order to align the appropriate slit
    with the neutron beam. These angles are not used in reduction.

    If the aperture is computed, then the value is stored
    in log key "source-aperture-diameter", with mili meter units

    Parameters
    ----------
    run: Mantid Run instance, MatrixWorkspace, file name, run number
    unit: str
        Length unit, either 'm' or 'mm'

    Returns
    -------
    float
        Source aperture diameter, in requested units
    """
    log_key = 'source-aperture-diameter'
    sl = SampleLogs(run)
    if log_key in sl.keys():
        sad = sl.single_value(log_key)  # units are 'mm'
    else:
        sad = source_aperture(run, unit='mm').diameter
        sl.insert(log_key, sad, unit='mm')
    if unit == 'm':
        sad /= 1000.0

    return sad


def insert_aperture_logs(ws):
    r"""
    Insert source and sample aperture diameters in the logs, as well as
    the distance between the source aperture and the sample. Units are in mm

    Parameters
    ----------
    ws: MatrixWorkspace
        Insert metadata in this workspace's logs
    """
    sl = SampleLogs(ws)
    if 'sample-aperture-diameter' not in sl.keys():
        sample_aperture_diameter(ws, unit='mm')  # function will insert the log
    if 'source-aperture-diameter' not in sl.keys():
        sad = source_aperture(ws, unit='mm').diameter
        sl.insert('source-aperture-diameter', sad, unit='mm')
    if 'source-aperture-sample-distance' not in sl.keys():
        sds = source_aperture(ws, unit='mm').distance_to_sample
        sl.insert('source-aperture-sample-distance', sds, unit='mm')
