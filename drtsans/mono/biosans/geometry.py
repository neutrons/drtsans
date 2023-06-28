# local imports
from drtsans.geometry import get_position_south_detector
from drtsans.mono.biosans.beam_finder import _calculate_neutron_drop
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.api import mtd
from mantid.api import Workspace as MantidWorkspace
from mantid.simpleapi import MoveInstrumentComponent, RotateInstrumentComponent
import numpy as np

# standard imports
from functools import lru_cache
from typing import Union

_reference_tubes = dict(
    wing_detector="bank49/tube4",  # wing tube closest to the beam
    midrange_detector="bank104/tube4",  # midrange tube closest to the beam
    detector1="bank1/tube1",
)  # south detector tube farthest from the beam

WING_RADIUS = 1.1633  # meters, radius of curvature of the wing detector
MIDRANGE_RADIUS = 4.0000  # meters, radius of curvature of the midrange detector
PHI_SPAN_MIDRANGE = 4.975  # degrees, horizontal angle (on the XZ plane) span of the midrange detector
#
TUBE_LENGTH = 1.046  # meters
REFERENCE_TUBES = dict(wing_detector="bank49/tube4", midrange_detector="bank104/tube4", detector1="bank1/tube1")
#
PIXELS_IN_TUBE = 256
PIXEL_HEIGHT = TUBE_LENGTH / PIXELS_IN_TUBE  # meters


def has_midrange_detector(input_workpace: Union[str, MantidWorkspace]) -> bool:
    r"""
    Check if the workspace has a midrange detector.

    Parameters
    ----------
    input_workpace
        Mantid workspace or name of the workspace

    Returns
    -------
    True if ``input_workspace`` has a midrange detector, False otherwise.
    """
    workspace = mtd[str(input_workpace)]
    return workspace.getInstrument().getComponentByName("midrange_detector") is not None


def set_position_south_detector(input_workspace, distance):
    r"""
    Position the south detector away from the origin of coordinates by some distance along the beam.

    This will also update log entry "sample_detector_distance" with the new distance of the detector from the sample.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object
    distance : float
        Distance from the sample to the south detector, in meters

    Raises
    ------
    ValueError : if distance is not greater than the distance from the wing detector to the sample, or
    the distance from the midrange detector to the sample if the midrange detector is present.
    """
    workspace = mtd[str(input_workspace)]
    mantid_instrument = workspace.getInstrument()
    midrange = mantid_instrument.getComponentByName(_reference_tubes["midrange_detector"])
    if midrange is not None:
        if distance <= MIDRANGE_RADIUS:
            raise ValueError(f"Invalid distance: {distance}. Must be greater than {MIDRANGE_RADIUS} meters.")
    else:
        if distance <= WING_RADIUS:
            raise ValueError(f"Invalid distance: {distance}. Must be greater than {WING_RADIUS} meters.")
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=distance, RelativePosition=False)
    sample = mantid_instrument.getSample()
    SampleLogs(workspace).insert("sample_detector_distance", distance - sample.getPos().Z())


def _get_angle_curved_detector(input_workspace, detector_name, offset_rotation, counter_clockwise=True):
    instrument = mtd[str(input_workspace)].getInstrument()
    curved_detector = instrument.getComponentByName(detector_name)
    angle_sign_flipper = 1.0 if counter_clockwise is True else -1.0
    return angle_sign_flipper * (curved_detector.getRotation().getEulerAngles()[0] - offset_rotation)


def get_angle_wing_detector(input_workspace):
    r"""
    Get the angle (in degrees) of the wing detector away from the beam axis. Increasing angles move
    the wing detector away from the beam axis.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object

    Returns
    -------
    float

    """
    OFFSET_ROTATION = -22.18  # (degrees) rotation of the wing detector such that it grazes the beam axis
    return _get_angle_curved_detector(input_workspace, "wing_detector", OFFSET_ROTATION, counter_clockwise=False)


def get_angle_midrange_detector(input_workspace):
    r"""
    Get the angle (in degrees) of the wing detector away from the beam axis. Increasing angles move
    the midrange detector away from the beam axis.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object

    Returns
    -------
    float

    """
    OFFSET_ROTATION = 2.48  # (degrees) rotation of the wing detector such that it grazes the beam axis
    return _get_angle_curved_detector(input_workspace, "midrange_detector", OFFSET_ROTATION, counter_clockwise=True)


def get_angle_south_detector(input_workspace):
    r"""
    Find the angle (in degrees) between the beam axis and the line joining the origin of coordinates and the edge
    of the south detector intersecting the horizontal (XZ) plane.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object

    Returns
    -------
    float
    """
    workspace = mtd[str(input_workspace)]
    half_span = workspace.getInstrument().getComponentByName("bank1/tube1").getPos().X()
    position = workspace.getInstrument().getComponentByName("detector1").getPos().Z()
    return float(np.degrees(np.arctan(half_span / position)))


def _set_angle_curved_panel(input_workspace, angle, detector_name, offset_rotation, counter_clockwise=True):
    if angle < 0.0:
        raise ValueError(f"Invalid angle: {angle}. Must be greater than 0 degrees.")
    workspace = mtd[str(input_workspace)]
    angle_sign_flipper = 1.0 if counter_clockwise is True else -1.0
    RotateInstrumentComponent(
        Workspace=workspace,
        ComponentName=detector_name,
        X=0,
        Y=1,
        Z=0,
        Angle=angle_sign_flipper * angle + offset_rotation,
        RelativeRotation=False,
    )


def set_angle_wing_detector(input_workspace, angle):
    r"""
    Position the wing detector at a given angle away from the beam axis.

    This will also update log entry "ww_rot_Readback" with the new angle.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object
    angle : float
        Angle from the beam axis, in degrees. Increasing angles move the wing detector away from the beam axis.

    Raises
    ------
    ValueError : if angle is not between 0 and 90 degrees
    """
    OFFSET_ROTATION = -22.18  # (degrees) rotation of the wing detector such that it grazes the beam axis
    SampleLogs(input_workspace).insert("ww_rot_Readback", angle)
    return _set_angle_curved_panel(input_workspace, angle, "wing_detector", OFFSET_ROTATION, counter_clockwise=False)


def set_angle_midrange_detector(input_workspace, angle):
    r"""
    Position the midrange detector at a given angle away from the beam axis.

    This will also update log entry "md_rot_Readback" with the new angle.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object
    angle : float
        Angle from the beam axis, in degrees. Increasing angles move the wing detector away from the beam axis.

    Raises
    ------
    ValueError : if angle is not between 0 and 90 degrees
    """
    OFFSET_ROTATION = 2.48  # (degrees) rotation of the wing detector such that it grazes the beam axis
    SampleLogs(input_workspace).insert("md_rot_Readback", angle)
    return _set_angle_curved_panel(
        input_workspace, angle, "midrange_detector", OFFSET_ROTATION, counter_clockwise=True
    )


def adjust_midrange_detector(input_workspace, criterium="fair_tube_shadowing"):
    r"""
    Find the optimal rotation angle on the horizontal (XZ plane) for the midrange detector according to some criterium.

    Criterium "fair_tube_shadowing" : the number of tubes in the south detector shadowed by the midrange detector
    coincides with the number of tubes in the mirrored-wing detector shadowing the midrange detector. The mirrored-wing
    detector is mirror image of the wing detector across the vertical (YZ) plane. Hence, if I(Q)_wing and I(Q)_midrange
    overlaps in N Q-values, I(Q)_south and I(Q)_midrange will also overlap in about N Q-values.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object
    criterium : str
        Criterium to use to find the optimal rotation angle. Valid values are:
        "fair_tube_shadowing"

    Returns
    -------
    float
    """
    # validate the input criterium
    if criterium != "fair_tube_shadowing":
        raise NotImplementedError(f"Invalid criterium: {criterium}")
    if criterium == "fair_tube_shadowing":
        phi_wing = get_angle_wing_detector(input_workspace)
        phi_south = get_angle_south_detector(input_workspace)
        # no fair shadowing is possible when the wing detector is rotated beyond:
        if phi_wing > phi_south + PHI_SPAN_MIDRANGE:
            phi = phi_south
        else:
            south_radius = get_position_south_detector(input_workspace)
            # the arc length of the wing detector shadowing the midrange detector when the midrange detector
            # is rotated by phi:
            #    arc_wing = (phi + phi_span_midrange - phi_wing) * radius_wing
            # arc length of the south detector shadowed by the phi-rotated midrange detector:
            #    arc_south = (phi_south - phi) * radius_south
            # equating arc_wing = arc_south and solving for phi gives:
            phi = ((phi_wing - PHI_SPAN_MIDRANGE) * WING_RADIUS + phi_south * south_radius) / (
                WING_RADIUS + south_radius
            )
        set_angle_midrange_detector(input_workspace, phi)
        return phi


@lru_cache(maxsize=128)
def midrange_to_wing_tubepixel(input_workspace):
    r"""
    Find the vertical position within a tube of the wing detector associated to a vertical position within a tube of
    the midrange detector, in pixel coordinates.

    This function takes into account the gravitational drop of the neutrons, and the fact that the midrange detector
    is positioned farther than the sample than the wing detector so that tubes in the midrange detector
    subtend a smaller solid angle than tubes in the wing detector.

    Typically, 70 pixels in the wing tube span the same vertical angle as the 256 pixels of the midrange tube.

    Parameters
    ----------
    input_workspace : str, ~mantid.simpleapi.Workspace
        Mantid workspace object

    Returns
    -------
    list
        List of integers, one for each pixel in the midrange detector, with the corresponding pixel
        in the wing detector
    """
    ws = mtd[str(input_workspace)]
    wavelength = SampleLogs(ws).single_value("wavelength")
    # valid approximation: use the origin of coordinates instead of the sample position
    wing_gravity_drop = _calculate_neutron_drop(WING_RADIUS, wavelength)
    midrange_gravity_drop = _calculate_neutron_drop(MIDRANGE_RADIUS, wavelength)
    #
    md2ww_pixel = list()
    l2_ratio = WING_RADIUS / MIDRANGE_RADIUS
    for midrange_pixel in range(256):
        midrange_y_coord = -TUBE_LENGTH / 2 - midrange_gravity_drop + midrange_pixel * PIXEL_HEIGHT
        wing_y_coord = l2_ratio * midrange_y_coord
        wing_pixel = int((wing_y_coord + TUBE_LENGTH / 2 + wing_gravity_drop) / PIXEL_HEIGHT)
        md2ww_pixel.append(wing_pixel)
    return md2ww_pixel


# quick id's for different components
info_ids = {
    "detector1": {"spectrum_info_range": (2, 49154)},  # 49154 is the first pixel in the next component
    "wing_detector": {"spectrum_info_range": (49154, 90114)},  # 90114 is the first pixel in the next component
    "midrange_detector": {"spectrum_info_range": (90114, 106498)},
}
