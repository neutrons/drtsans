from scipy import constants

import drtsans.beam_finder as bf
from mantid import mtd
from mantid.kernel import logger
from drtsans.samplelogs import SampleLogs


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
    sdd_main_detector = ws.getInstrument().getComponentByName('detector1').getPos().Z()
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


def find_beam_center(input_workspace, method='center_of_mass', mask=None,
                     sdd_wing_detector=1.13, **kwargs):
    """Finds the beam center in a 2D SANS data set.
    Developer: Ricardo Ferraz Leal <rhf@ornl.gov>

    based on drtsans.beam_finder.find_beam_center

    Parameters
    ----------
    input_workspace: str, Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str or list
        Tubes to be masked
    sdd_wing_detector : float, optional
        Sample Detector Distance,  in meters, of the wing detector,
        by default 1.13
    kwargs: dict
        Parameters to be passed to the method to calculate the center

    Returns
    -------
    Tuple of 3 floats
        center_x, center_y, center_y corrected for gravity
        center_y it is usually used to correct BIOSANS wing detector
        Y position.
    """
    center_x, center_y = bf.find_beam_center(input_workspace, method, mask, **kwargs)

    center_y_gravity = _beam_center_gravitational_drop(
        input_workspace, center_y, sdd_wing_detector)

    logger.information("Beam Center: x={:.3} y={:.3} y_gravity={:.3}.".format(
        center_x, center_y, center_y_gravity))
    return center_x, center_y, center_y_gravity


def center_detector(input_workspace, center_x, center_y, center_y_gravity):
    """Center the detector and adjusts the heigh for the wing detector

    **Mantid algorithms used:**
    :ref:`MoveInstrumentComponent <algm-MoveInstrumentComponent-v1>`,

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
    # move the main detector
    bf.center_detector(input_workspace, center_x, center_y)

    # movethe wing detector for the gravity drop
    bf.center_detector(input_workspace, center_x, center_y_gravity, component='wing_detector')
