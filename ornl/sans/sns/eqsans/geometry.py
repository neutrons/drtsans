from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.simpleapi import MoveInstrumentComponent

from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (sample_source_distance,
                                sample_detector_distance,
                                detector_name)

__all__ = ['detector_z_log',
           'translate_sample_by_z', 'translate_detector_by_z',
           'center_detector', 'sample_aperture_diameter',
           'source_aperture_diameter']

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
    MoveInstrumentComponent(ws, Z=z,
                            ComponentName='sample-position',
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
        sl = SampleLogs(ws)
        z = 1e-3 * sl.single_value(detector_z_log)  # assumed in mili-meters

    kwargs = dict(ComponentName=detector_name(ws),
                  RelativePosition=relative)
    MoveInstrumentComponent(ws, Z=z, **kwargs)


def translate_detector_by_z(ws, z, **kwargs):
    r"""
    Simplified call to translate_detector_z
    """
    return translate_detector_z(ws, z=z, relative=True, **kwargs)


def center_detector(ws, x, y, units='m', relative=False):
    r"""
    Move the detector on the XY plane.

    Usually `x` and `y` will be the absolute coordinates of the beam
    impinging on the detector.

    Parameters
    ----------
    ws: Workspace
        Input workspace containing the instrument
    x: float
        Final position or translation along the X-axis
    y: float
        Final position or translation along the Y-axis
    units: str
        Either meters 'm' or mili-meters 'mm'
    relative: Bool
        Apply translation if True, otherwise `x` and `y` are final coordinates

    Returns
    =======
    numpy.ndarray
        Detector vector position
    """
    t_x = x if units == 'm' else x / 1.e3
    t_y = y if units == 'm' else y / 1.e3
    t_z = 0.0
    if relative is False:
        i = ws.getInstrument()
        # preserve the Z coordinate value
        t_z = i.getComponentByName(detector_name(i)).getPos()[-1]
    MoveInstrumentComponent(ws, X=t_x, Y=t_y, Z=t_z,
                            ComponentName=detector_name(ws),
                            RelativePosition=relative)

    # Recalculate distance from sample to detector
    sdd = sample_detector_distance(ws, units='mm', search_logs=False)
    SampleLogs(ws).insert('sample-detector-distance', sdd, unit='mm')
    instrument = ws.getInstrument()
    det = instrument.getComponentByName(detector_name(instrument))
    return np.array(det.getPos())


def sample_aperture_diameter(run, unit='mm'):
    r"""
    Find the sample aperture diameter from the logs.

    Log keys searched are 'sample-aperture-diameter' and 'beamslit4'.

    Parameters
    ----------
    run: Mantid Run instance, MatrixWorkspace, file name, run number
        Input from which to find the aperture
    unit: str
        return aperture in requested length unit, either 'm' or 'mm'

    Returns
    -------
    float
        Sample aperture diameter, in requested units
    """
    sl = SampleLogs(run)
    for log_key in ('sample-aperture-diameter', 'beamslit4'):
        if log_key in sl.keys():
            sad = sl.single_value(log_key)
            break
    if 'sample-aperture-diameter' not in sl.keys():
        sl.insert('sample-aperture-diameter', sad, unit='mm')
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
    sl = SampleLogs(other)
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
        sample_aperture_diameter(ws, unit='mm')  # this function will insert the log
    if 'source-aperture-diameter' not in sl.keys():
        sad = source_aperture(ws, unit='mm').diameter
        sl.insert('source-aperture-diameter', sad, unit='mm')
    if 'source-aperture-sample-distance' not in sl.keys():
        sds = source_aperture(ws, unit='mm').distance_to_sample
        sl.insert('source-aperture-sample-distance', sds, unit='mm')
