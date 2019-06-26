from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import (MoveInstrumentComponent, AddSampleLog)

from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogsReader
from ornl.sans.geometry import sample_source_distance

__all__ = ['detector_name', 'detector_z_log',
           'translate_sample_by_z', 'translate_detector_by_z']

detector_name = 'detector1'
detector_z_log = 'detectorZ'


def translate_sample_by_z(ws, z):
    r"""
    Shift the position of the sample by the desired amount

    Parameters
    ----------
    ws: Workspace
        Input workspace containing instrument file
    z: float
        Translation to be applied
    """
    MoveInstrumentComponent(ws, Z=z, ComponentName='sample',
                            RelativePosition=True)


def translate_detector_z(ws, z=None, relative=True):
    r"""
    Adjust Z-coordinate of detector bank in instrument file.


    Parameters
    ----------
    ws: Workspace
        Input workspace containing instrument file
    z: float
        Translation to be applied, in units of meters. If `None`, log_key
        stored in `detector_z_log` is used
    relative: bool
        If True, add to the current z-coordinate. If False, substitute
        the current z-coordinate with the new value.
    """
    if z is None:
        sl = SampleLogsReader(ws)
        z = 1e-3 * sl.single_value(detector_z_log)  # assumed in mili-meters

    kwargs = dict(ComponentName=detector_name,
                  RelativePosition=relative)
    MoveInstrumentComponent(ws, Z=z, **kwargs)


def translate_detector_by_z(ws, z, **kwargs):
    r"""
    Simplified call to translate_detector_z
    """
    return translate_detector_z(ws, z=z, relative=True, **kwargs)


def sample_aperture_diameter(other, unit='mm'):
    r"""
    Find the sample aperture diameter

    Parameters
    ----------
    other: Run, MatrixWorkspace, file name, run number
    unit: str
        Length unit, either 'm' or 'mm'
    Returns
    -------
    float
        Sample aperture diameter, in requested units
    """
    sl = SampleLogsReader(other)
    sad = sl.single_value('beamslit4')
    if unit == 'm':
        sad /= 1000.0
    return sad


@namedtuplefy
def source_aperture(other, unit='mm'):
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
    sl = SampleLogsReader(other)
    slit_indexes = [int(sl[log_key].value.mean()) - 1 for log_key in
                    ['vBeamSlit', 'vBeamSlit2', 'vBeamSlit3']]
    diameter = 20.0  # default slit size
    asd = -1.0  # aperture to sample distance
    ssd = sample_source_distance(other, units='mm')
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


def source_aperture_diameter(other, unit='mm'):
    r"""
    Find the source aperture diameter

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
    float
        Source aperture diameter, in requested units
    """
    return source_aperture(other, unit=unit).diameter


def insert_aperture_logs(ws):
    r"""
    Insert source and sample aperture diameters in the logs.
    Units are in mm

    Parameters
    ----------
    ws: MatrixWorkspace
        Insert metadata in this workspace's logs
    """
    kw = dict(LogType='Number', LogUnit='mm', NumberType='Double')
    AddSampleLog(ws, Logname='sample-aperture-diameter',
                 LogText=str(sample_aperture_diameter(ws)), **kw)
    sa = source_aperture(ws)
    AddSampleLog(ws, Logname='source-aperture-diameter',
                 LogText=str(sa.diameter), **kw)
    AddSampleLog(ws, Logname='source-aperture-sample-distance',
                 LogText=str(sa.distance_to_sample), **kw)

