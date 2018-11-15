from mantid.simpleapi import (
    SANSMaskDTP, FindCenterOfMassPosition)
from mantid.kernel import logger


def _beam_center_gravitational_drop(ws, beam_center_y, sdd=1.13):
    '''
    This method is used for correcting for gravitational drop
    ws - it's a workspace where the
    sdd - sample detector distance to apply the beam center
    '''
    def calculate_neutron_drop(path_length, wavelength):
        '''
        Calculate the gravitational drop of the neutrons
        path_length in meters
        wavelength in Angstrom
        Return the Y drop in meters
        '''
        wavelength *= 1e-10
        neutron_mass = 1.674927211e-27
        gravity = 9.80665
        h_planck = 6.62606896e-34
        l_2 = (gravity * neutron_mass**2 /
               (2.0 * h_planck**2)) * path_length**2
        return wavelength**2 * l_2

    logger.information("Beam Center Y before: %.2f meters" % beam_center_y)

    i = ws.getInstrument()
    # distance from the sample to the main detector
    distance_detector1 = i.getComponentByName("detector1").getPos()[2]
    path_length = distance_detector1 - sdd
    logger.debug("SDD detector1 = %.3f meters. SDD for wing = %.3f meters." % (
        distance_detector1, sdd))
    logger.debug("Path length for gravitational drop = %.3f meters." %
                 (path_length))
    r = ws.run()
    wavelength = r.getProperty("wavelength").value
    logger.debug("Wavelength = %.2f A." % (wavelength))
    drop = calculate_neutron_drop(path_length, wavelength)
    logger.debug("Gravitational drop = %.6f meters." % (drop))
    new_beam_center_y = beam_center_y + drop
    logger.information("Beam Center Y after: %.2f meters" % new_beam_center_y)
    return new_beam_center_y


def direct_beam_center(input_ws, tubes_to_mask=None, sdd_wing_detector=1.13):
    '''
    Return beam center x, y in meters
    '''

    if tubes_to_mask is not None:
        SANSMaskDTP(InputWorkspace=input_ws, Tube=tubes_to_mask)

    center = FindCenterOfMassPosition(InputWorkspace=input_ws)
    center_x, center_y = center

    center_y_gravity = _beam_center_gravitational_drop(
        input_ws, center_y, sdd_wing_detector)

    return center_x, center_y, center_y_gravity
