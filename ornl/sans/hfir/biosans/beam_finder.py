from mantid.simpleapi import FindCenterOfMassPosition
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


def direct_beam_center(input_workspace, center_x=0, center_y=0,
                       tolerance=0.00125, direct_beam=True, beam_radius=0.0155,
                       sdd_wing_detector=1.13):
    """Finds the beam center in a 2D SANS data set.

    Parameters
    ----------
    input_workspace : MatrixWorkspace, str
        The beamcenter workspace
    center_x : int, optional
        Estimate for the X beam center in meters, by default 0
    center_y : int, optional
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
    Tuple of 3 ints
        center_x, center_y, center_y corrected for gravity
        center_y it is usually used to correct BIOSANS wing detector
        Y position.
    """

    center = FindCenterOfMassPosition(
        InputWorkspace=input_workspace, CenterX=center_x, CenterY=center_y,
        Tolerance=tolerance, DirectBeam=direct_beam, BeamRadius=beam_radius)

    center_x, center_y = center

    center_y_gravity = _beam_center_gravitational_drop(
        input_workspace, center_y, sdd_wing_detector)

    logger.information("Beam Center: x={:.3} y={:.3} y_gravity={:.3}".format(
        center_x, center_y, center_y_gravity))
    return center_x, center_y, center_y_gravity
