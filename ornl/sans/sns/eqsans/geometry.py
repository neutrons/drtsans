from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.simpleapi import (MoveInstrumentComponent, Integration,
                              FindCenterOfMassPosition, LoadMask,
                              MaskDetectors)

from ornl.settings import (namedtuplefy,
                           unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (sample_source_distance,
                                sample_detector_distance,
                                detector_name)

__all__ = ['detector_z_log',
           'translate_sample_by_z', 'translate_detector_by_z',
           'center_detector', 'direct_beam_center', 'sample_aperture_diameter',
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
    MoveInstrumentComponent(Workspace=str(ws), Z=z,
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

    MoveInstrumentComponent(Workspace=ws, Z=z, ComponentName=detector_name(ws),
                            RelativePosition=relative)


def translate_detector_by_z(ws, z, **kwargs):
    r"""
    Simplified call to translate_detector_z
    """
    return translate_detector_z(ws, z=z, relative=True, **kwargs)


def center_detector(ws, method='center_of_mass',
                    x=None, y=None, units='m', relative=False,
                    move_detector=True):
    r"""
    Move the detector on the XY plane to center the beam location

    An estimation of the center of the detector is carried out
    according to `method`, unless `x` and `y` absolute coordinates
    are provided (`relative=False`) in which case `x` and `y` are used.
    If `x` and `y` are a translation (`relative=False`) then both
    `method` and the translation will be applied.

    case 1: position detector after a center of mass estimation
        center_detector(ws)
    case 2: position detector at absolute coordinates (x0,y0)
        center_detector(ws, x=x0, y=y0)
    case 3: translate the detector by (x0, y0)
        enter_detector(ws, method=None, x=x0, y=y0, relative=True)
    case 4: position detector after a center of mass estimation followed
            by a translation (x0, y0)
         center_detector(ws, x=x0, y=y0, relative=True)

    Parameters
    ----------
    ws: Workspace
        Input workspace containing the instrument
    method: str
        Method to estimate the center of the beam. `None` for no method
    x: float
        Final position or translation along the X-axis
    y: float
        Final position or translation along the Y-axis
    units: str
        units of `x` and `y`. Either meters 'm' or mili-meters 'mm'
    relative: Bool
        Values of `x` and `y` are either absolute coordinates or a
        translation.
    move_detector: bool
        Only calculate the final position if this is False

    Returns
    =======
    numpy.ndarray
        Detector vector position
    """
    i = ws.getInstrument()
    rs = i.getComponentByName(detector_name(i)).getPos()
    rf = np.zeros(3)

    # Use `method` when we are not passing absolute coordinates
    abs_xy = x is not None and y is not None and relative is False
    if method is not None and abs_xy is False:
        method_to_alg = dict(center_of_mass=FindCenterOfMassPosition)
        ws_flattened = Integration(InputWorkspace=ws, OutputWorkspace=uwd())
        t_x, t_y = list(method_to_alg[method](InputWorkspace=ws_flattened))
        rs = np.array([t_x, t_y, rs[-1]])
        rf = rs

    if x is not None and y is not None:
        t_x = x if units == 'm' else x / 1.e3
        t_y = y if units == 'm' else y / 1.e3
        if relative is True:
            rf = rs + np.array([t_x, t_y, 0.0])
        else:
            rf = np.array([t_x, t_y, rs[-1]])

    # Recalculate distance from sample to detector
    if move_detector is True:
        MoveInstrumentComponent(ws, X=rf[0], Y=rf[1], Z=rf[2],
                                ComponentName=detector_name(ws),
                                RelativePosition=False)
        sdd = sample_detector_distance(ws, units='mm', search_logs=False)
        SampleLogs(ws).insert('sample-detector-distance', sdd, unit='mm')

    return rf


def find_beam_center(input_workspace, method='center_of_mass', mask=None):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    Parameters
    ----------
    input_workspace: Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str, MaskWorkspace
        Path to mask file, or MaskWorkspace object

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters)
    """
    if isinstance(mask, str):
        mask = LoadMask(mask, OutputWorkspace=uwd())
    if mask is not None:
        w = MaskDetectors(input_workspace, OutputWorkspace=uwd())
    r = center_detector(input_workspace, method=method, move_detector=False)
    return (r[0], r[1])


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
        sample_aperture_diameter(ws, unit='mm')  # function will insert the log
    if 'source-aperture-diameter' not in sl.keys():
        sad = source_aperture(ws, unit='mm').diameter
        sl.insert('source-aperture-diameter', sad, unit='mm')
    if 'source-aperture-sample-distance' not in sl.keys():
        sds = source_aperture(ws, unit='mm').distance_to_sample
        sl.insert('source-aperture-sample-distance', sds, unit='mm')
