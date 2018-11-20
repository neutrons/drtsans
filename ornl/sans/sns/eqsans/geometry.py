from __future__ import (absolute_import, division, print_function)

from ornl.settings import namedtuplefy
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import sample_source_distance


def sample_aperture_diameter(other):
    r"""
    Find the sample aperture diameter

    Parameters
    ----------
    other: Run, MatrixWorkspace, file name, run number

    Returns
    -------
    float
        Sample aperture diameter, in mili-meters
    """
    sl = SampleLogs(other)
    return float(sl.beamslit4.value.mean())


@namedtuplefy
def source_aperture(other):
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

    Returns
    -------
    namedtuple
        Fields of the name tuple
        - float: diameter, in mili-meters
        - float: distance to sample, in mili-meters
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
    return dict(diameter=diameter, distance_to_sample=asd)


def source_aperture_diameter(other):
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

    Returns
    -------
    float
        Source aperture diameter, in mili-meters
    """
    return source_aperture(other).diameter


def insert_aperture_logs(ws):
    r"""
    Insert source and sample aperture diameters in the logs.

    Parameters
    ----------
    ws: MatrixWorkspace
        Insert metadata in this workspace's logs
    """
    sl = SampleLogs(ws)
    sl['sample_aperture-diameter'] = sample_aperture_diameter(ws)
    sa = source_aperture(ws)
    sl['source_aperture-diameter'] = sa.diameter
    sl['source_aperture-sample-distance'] = sa.distance_to_sample
