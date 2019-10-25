from scipy import constants

import drtsans.beam_finder as bf
from mantid import mtd
from mantid.kernel import logger
from drtsans.samplelogs import SampleLogs

__all__ = ['center_detector', 'find_beam_center']


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


def _beam_center_gravitational_drop(ws, beam_center_y, sample_det_cent_wing_detector=1.13):
    """This method is used for correcting for gravitational drop.
    It is used in the biosans for the correction of the beam center
    in the wing detector. The center in the wing detector will be higher
    because the neutrons fall due to gravity until they hit the main detector.

    Parameters
    ----------
    ws : ~mantid.api.MatrixWorkspace, str
        The workspace to get the sample to the main detector distance
    beam_center_y : float
        The Y center coordinate in the main detector in meters
    sample_det_cent_wing_detector : float, optional
        :ref:`sample to detector center distance <devdocs-standardnames>` of the wing detector
        in meters to apply the beam center

    Returns
    -------
    float
        The new y beam center corrected for distance
    """

    ws = mtd[str(ws)]

    sl = SampleLogs(ws)
    sample_det_cent_main_detector = ws.getInstrument().getComponentByName('detector1').getPos().Z()
    path_length = sample_det_cent_main_detector - sample_det_cent_wing_detector

    wavelength = sl.wavelength.value

    r = ws.run()
    wavelength = r.getProperty("wavelength").value

    drop = _calculate_neutron_drop(path_length, wavelength)
    new_beam_center_y = beam_center_y - drop
    logger.information("Beam Center Y before gravity (drop = {:.3}): {:.3}"
                       " after = {:.3}.".format(
                            drop, beam_center_y, new_beam_center_y))

    return new_beam_center_y


def find_beam_center(input_workspace, method='center_of_mass', mask=None, mask_options={}, centering_options={},
                     sample_det_cent_wing_detector=1.13):
    """Finds the beam center in a 2D SANS data set.
    This is based on (and uses) :func:`drtsans.find_beam_center`

    Developer: Ricardo Ferraz Leal <rhf@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
    method: str
        Method to calculate the beam center. Available methods are:
        - 'center_of_mass', invokes :ref:`FindCenterOfMassPosition <algm-FindCenterOfMassPosition-v1>`.
    mask: mask file path, `MaskWorkspace``, :py:obj:`list`.
        Mask to be passed on to ~drtsans.mask_utils.mask_apply.
    mask_options: dict
        Additional arguments to be passed on to ~drtsans.mask_utils.mask_apply.
    centering_options: dict
        Arguments to be passed on to the centering method.
    sample_det_cent_wing_detector : float
        :ref:`sample to detector center distance <devdocs-standardnames>`,
        in meters, of the wing detector.

    Returns
    -------
    tuple
        Three float numbers:
        ``(center_x, center_y, center_y)`` corrected for gravity.
        ``center_y`` is usually used to correct BIOSANS wing detector Y position.
    """
    center_x, center_y = bf.find_beam_center(input_workspace, method, mask,
                                             mask_options=mask_options, centering_options=centering_options)

    center_y_gravity = _beam_center_gravitational_drop(
        input_workspace, center_y, sample_det_cent_wing_detector)

    logger.information("Beam Center: x={:.3} y={:.3} y_gravity={:.3}.".format(
        center_x, center_y, center_y_gravity))
    return center_x, center_y, center_y_gravity


def center_detector(input_workspace, center_x, center_y, center_y_gravity):
    """Center the detector and adjusts the height for the wing detector
    This uses :func:`drtsans.center_detector`

    **Mantid algorithms used:**
    :ref:`MoveInstrumentComponent <algm-MoveInstrumentComponent-v1>`,

    Parameters
    ----------
    input_workspace : ~mantid.api.MatrixWorkspace, str
        The workspace to be centered
    center_x : float
        in meters
    center_y : float
        in meters
    center_y_gravity : float
        in meters

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        reference to the corrected ``input_workspace``
    """
    # move the main detector
    bf.center_detector(input_workspace, center_x, center_y)

    # movethe wing detector for the gravity drop
    bf.center_detector(input_workspace, center_x, center_y_gravity, component='wing_detector')
