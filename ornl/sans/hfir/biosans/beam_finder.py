from scipy import constants

from mantid import mtd
from mantid.kernel import logger
# https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
# https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html
from mantid.simpleapi import FindCenterOfMassPosition, MoveInstrumentComponent
from ornl.sans.samplelogs import SampleLogs


def _calculate_neutron_drop(path_length, wavelength):
    """ Calculate the gravitational drop of the neutrons

    Parameters
    ----------
    path_length : float
        path_length in meters
    wavelength : float
        wavelength in Angstrom

    Returns
    -------
    float
        Return the Y drop in meters
    """
    wavelength *= 1e-10
    neutron_mass = constants.neutron_mass
    gravity = constants.g
    h_planck = constants.Planck
    l_2 = (gravity * neutron_mass**2 /
           (2.0 * h_planck**2)) * path_length**2
    return wavelength**2 * l_2


def _beam_center_gravitational_drop(ws, beam_center_y, sdd_wing_detector=1.13):
    """This method is used for correcting for gravitational drop.
    It is used in the biosans for the correction of the beam center
    in the wing detector. The center in the wing detector will be higher
    because the neutrons fall due to gravity until they hit the main detector.

    Parameters
    ----------
    ws : Workspace2D, str
        The workspace to get the sample to the main detector distance
    beam_center_y : float
        The Y center coordinate in the main detector in meters
    sdd : float, optional
        sample detector distance in meters to apply the beam center,
        by default 1.13

    Returns
    -------
    float
        The new y beam center corrected for distance
    """

    ws = mtd[str(ws)]

    sl = SampleLogs(ws)
    sdd_main_detector = sl['sample-detector-distance'].value*1e-3  # meters
    path_length = sdd_main_detector - sdd_wing_detector

    wavelength = sl.wavelength.value

    r = ws.run()
    wavelength = r.getProperty("wavelength").value

    drop = _calculate_neutron_drop(path_length, wavelength)
    new_beam_center_y = beam_center_y - drop
    logger.information("Beam Center Y before gravity (drop = {:.3}): {:.3}"
                       " after = {:.3}.".format(
                            drop, beam_center_y, new_beam_center_y))

    return new_beam_center_y


def direct_beam_center(input_workspace, center_x_estimate=0,
                       center_y_estimate=0, tolerance=0.00125,
                       direct_beam=True, beam_radius=0.0155,
                       sdd_wing_detector=1.13):
    """Finds the beam center in a 2D SANS data set.
    Developer: Ricardo Ferraz Leal <rhf@ornl.gov>

    Parameters
    ----------
    input_workspace : MatrixWorkspace, str
        The beamcenter workspace
    center_x_estimate : float, optional
        Estimate for the X beam center in meters, by default 0
    center_y_estimate : float, optional
        Estimate for the Y beam center in meters, by default 0
    tolerance : float, optional
        Tolerance on the center of mass position between each iteration in m,
        by default 0.00125
    direct_beam : bool, optional
        If true, a direct beam calculation will be performed. Otherwise, the
        center of mass of the scattering data will be computed by excluding
        the beam area., by default True
    beam_radius : float, optional
        Radius of the beam area, in meters, used the exclude the beam when
        calculating the center of mass of the scattering pattern,
        by default 0.0155
    sdd_wing_detector : float, optional
        Sample Detector Distance,  in meters, of the wing detector,
        by default 1.13

    Returns
    -------
    Tuple of 3 floats
        center_x, center_y, center_y corrected for gravity
        center_y it is usually used to correct BIOSANS wing detector
        Y position.
    """

    center = FindCenterOfMassPosition(
        InputWorkspace=input_workspace, CenterX=center_x_estimate,
        CenterY=center_y_estimate, Tolerance=tolerance, DirectBeam=direct_beam,
        BeamRadius=beam_radius)

    center_x, center_y = center

    center_y_gravity = _beam_center_gravitational_drop(
        input_workspace, center_y, sdd_wing_detector)

    logger.information("Beam Center: x={:.3} y={:.3} y_gravity={:.3}.".format(
        center_x, center_y, center_y_gravity))
    return center_x, center_y, center_y_gravity


def beam_center(center_x, center_y, wavelength, sdd_main_detector,
                sdd_wing_detector=1.13):
    """Given the coordinates of the beam center calculate the beam center
    of the Y coordinate in the wing detector.

    Parameters
    ----------
    center_x : float
        beam center X in meters found in the main detector
    center_y : float
        beam center y in meters found in the main detector
    wavelength : float
        in Angstroms
    sdd_main_detector : float
        in meters
    sdd_wing_detector : float, optional
        in meters, by default 1.13
    """

    drop = _calculate_neutron_drop(sdd_main_detector - sdd_wing_detector,
                                   wavelength)
    center_y_gravity = center_y + drop

    logger.information("Beam Center: x={:.3} y={:.3} y_gravity={:.3}".format(
        center_x, center_y, center_y_gravity))
    return center_x, center_y, center_y_gravity


# API

def center_detector(input_workspace, center_x, center_y, center_y_gravity):
    """Center the detector and adjusts the heigh for the wing detector

    Parameters
    ----------
    input_workspace : Workspace2D, str
        The workspace to be centered
    center_x : float
        in meters
    center_y : float
        in meters
    center_y_gravity : float
        in meters

    Returns
    -------
    Workspace2D
        reference to the input_workspace
    """

    MoveInstrumentComponent(
        Workspace=input_workspace, ComponentName='detector1',
        X=-center_x, Y=-center_y)

    # Now let's correct the wing detector for the gravity drop
    # Relative movement up words
    MoveInstrumentComponent(
        Workspace=input_workspace, ComponentName='wing_detector',
        X=-center_x,
        Y=-center_y_gravity)

    return input_workspace
