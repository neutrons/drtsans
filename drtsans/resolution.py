import numpy as np
from scipy import constants


__all__ = ['InstrumentSetupParameters', 'calculate_sigma_theta_prefactor', 'calculate_sigma_geometry',
           'calculate_sigma_theta_geometry', 'calculate_sigma_theta_gravity']


class InstrumentSetupParameters(object):
    """
    Class to contain the parameters used to calculate Q resolution
    """
    def __init__(self, l1, sample_det_center_dist, source_aperture_radius, sample_aperture_radius,
                 pixel_size_x, pixel_size_y):
        """
        Initialization to set all the parameters (6) to calculate momentrum transfer resolution

        Parameters
        ----------
        l1:
            source to sample distance
        sample_det_center_dist:
            sample detector (bank) center distance
        source_aperture_radius:
            source aperture radius (meter)
        sample_aperture_radius:
            sample aperture radius (meter)
        pixel_size_x:
            pixel linear size along X direction (meter)
        pixel_size_y:
            pixel linear size along Y direction (meter)
        """
        self._l1 = l1
        self._sample_det_center_dist = sample_det_center_dist
        self._source_aperture = source_aperture_radius
        self._sample_aperture = sample_aperture_radius
        self._pixel_size_x = pixel_size_x
        self._pixel_size_y = pixel_size_y

    def __str__(self):
        """
        Nice output string
        :return:
        """
        out = 'L1 = {} (m)\nSample-Detector-Center-Distance (L2)= {} (m)\n' \
              ''.format(self.l1, self._sample_det_center_dist)
        out += 'Source aperture radius (R1) = {} (m)\n'.format(self._source_aperture)
        out += 'Sample aperture radius (R2) = {} (m)\n'.format(self._sample_aperture)
        out += 'Pixel size = {}, {} (m, m)'.format(self._pixel_size_x, self._pixel_size_y)

        return out

    @property
    def l1(self):
        """
        Get L1 value
        :return: L1 (meter)
        """
        return self._l1

    @property
    def sample_det_center_distance(self):
        """
        Distance from sample to detector bank center,
        which is L2 in the SANS master document
        :return: sample detector center distance, aka SANS L2 (meter)
        """
        return self._sample_det_center_dist

    @property
    def source_aperture_radius(self):
        """
        Source aperture radius, which is R1 in SANS master document
        :return: source aperture radius (R1) in meter
        """
        return self._source_aperture

    @property
    def sample_aperture_radius(self):
        """
        Sample aperture radius, which is R2 in SANS master document
        :return: sample aperture radius (R2) in meter
        """
        return self._sample_aperture

    @property
    def pixel_size_x(self):
        """
        Detector pixel size along X direction
        :return: detector size along X direction in meter
        """
        return self._pixel_size_x

    @property
    def pixel_size_y(self):
        """
        Detector pixel size along Y direction
        :return: detector size along Y direction in meter
        """
        return self._pixel_size_y


def calculate_sigma_theta_prefactor(wavelength, pixel_info, instrument_parameters):
    r"""
    Calculates for every pixel and wavelength

    .. math::

       \left(\frac{2\pi\cos\theta\cos^2(2\theta)}{\lambda L_2}\right)^2


    Parameters
    ----------

    wavelength: ~np.array
        the array of wavelengths (same shape as momentum transfer)
    pixel_info: ~collections.namedtuple
        A namedtuple with fields for two_theta, azimuthal, l2, keep
    instrument_parameters: InstrumentSetupParameters
        Information abot instrument

    Returns
    -------
    float
        The coefficient described above
    """
    two_theta = pixel_info.two_theta.reshape(-1, 1)
    L2 = instrument_parameters.sample_det_center_distance
    return np.square(2 * np.pi * np.cos(0.5 * two_theta) * np.cos(two_theta)**2 / wavelength / L2)


def calculate_sigma_theta_geometry(instrument_parameters, mode):
    r"""
    Calculates

    .. math::

       \left(\frac {L_2}{L_1}\right)^2\frac{R_1^2}{4}+\left(\frac {L_1+L_2}{L_1}\right)^2\frac{R_2^2}{4}+
       \frac {1}{12}(\Delta R)^2

    If mode is "scalar", :math:`((\Delta R)^2=(\Delta x)^2+(\Delta y)^2)/2`, else

    :math:`(\Delta R)^2=[(\Delta x)^2,(\Delta y)^2]`. The formula for scalar is consistent with
    the equations 10.3 and 10.4 in the master document. when you add the two together, the geometry
    part is twice the contribution of :math:`(\Delta R)^2` plus the gravity part.


    Parameters
    ----------

    instrument_parameters: InstrumentSetupParameters
        Information abot instrument
    mode: str
        One of "scalar", "azimuthal", "crystalographic"

    Returns
    -------
    float or list
        The coefficient described above
    """
    L1 = instrument_parameters.l1
    L2 = instrument_parameters.sample_det_center_distance
    R1 = instrument_parameters.source_aperture_radius
    R2 = instrument_parameters.sample_aperture_radius
    dx2 = np.square(instrument_parameters.pixel_size_x)
    dy2 = np.square(instrument_parameters.pixel_size_y)

    if mode == "scalar":
        pixel_size2 = 0.5 * (dx2 + dy2)
    elif mode == "azimuthal":
        pixel_size2 = np.array([dx2, dy2])
    return 0.25 * np.square(L2 / L1 * R1)  \
        + 0.25 * np.square((L1 + L2) / L1 * R2) + pixel_size2 / 12.0


def calculate_sigma_theta_gravity(instrument_parameters, wavelength, delta_wavelength):
    r"""
    Calculates

    .. math::

       \frac 23 B^2\lambda^2(\Delta\lambda)^2

    where :math:`B=g m_N^2L_2(L_1+L_2)/(2h^2)`

    """
    # derived constant where:
    # h = 6.626e-34    # m^2 kg s^-1
    # m_n = 1.675e-27  # kg
    # g = 9.8          # m s^-2
    G_MN2_OVER_H2 = constants.g * np.square(constants.neutron_mass / constants.h)  # Unit as m, s, Kg
    L1 = instrument_parameters.l1
    L2 = instrument_parameters.sample_det_center_distance
    B = 0.5 * G_MN2_OVER_H2 * L2 * (L1 + L2) * 1.E-20
    return 2. * np.square(B * wavelength * delta_wavelength) / 3.


def calculate_sigma_geometry(mode, wavelength, delta_wavelength, pixel_info, instrument_parameters):
    r"""
    Calculates the Q independent part of the resolution, the common parts in formula 10.3 - 10.6

    Parameters
    ----------
    mode: str
        One of "scalar", "azimuthal", "crystalographic"
    wavelength: ~np.array
        the array of wavelengths (same shape as momentum transfer)
    delta_wavelength: ~np.array
        the array of wavelength widths (same shape as momentum transfer)
    pixel_info: ~collections.namedtuple
        A namedtuple with fields for two_theta, azimuthal, l2, keep
    instrument_parameters: InstrumentSetupParameters
        Information abot instrument


    """
    factor = calculate_sigma_theta_prefactor(wavelength, pixel_info, instrument_parameters)
    geometry_part = calculate_sigma_theta_geometry(instrument_parameters, mode)
    gravity_part = calculate_sigma_theta_gravity(instrument_parameters, wavelength, delta_wavelength)

    if mode == "scalar":
        return factor * (geometry_part * 2 + gravity_part)
    if mode == "azimuthal":
        return [factor * geometry_part[0], factor * (geometry_part[1] + gravity_part)]
