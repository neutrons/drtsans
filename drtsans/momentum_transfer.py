"""
Resolution calculations common to all SANS
"""
import numpy as np
from scipy import constants
import collections

# derived constant where:
# h = 6.626e-34    # m^2 kg s^-1
# m_n = 1.675e-27  # kg
# g = 9.8          # m s^-2
G_MN2_OVER_H2 = constants.g * np.square(constants.neutron_mass / constants.h)  # Unit as m, s, Kg

""" Named tuple for momentum transfer and Q resolution
"""
MomentumTransfer = collections.namedtuple('MomentumTransfer', 'q qx qy qz dqx dqy')


class InstrumentSetupParameters(object):
    """
    Class to contain the parameters used to calculate Q resolution
    """
    def __init__(self, l1, sample_det_center_dist, source_aperture_radius, sample_aperture_radius,
                 pixel_size_x, pixel_size_y):
        """
        Initialization to set all the parameters (6) to calculate momentrum transfer resolution
        :param l1: L1 (source to sample)
        :param sample_det_center_dist: sample detector (bank) center distance
        :param source_aperture_radius: source aperture radius (meter)
        :param sample_aperture_radius: sample aperture radius (meter)
        :param pixel_size_x: pixel linear size along X direction (meter)
        :param pixel_size_y: pixel linear size along Y direction (meter)
        """
        self._l1 = l1
        self._sample_det_center_dist = sample_det_center_dist
        self._source_aperture = source_aperture_radius
        self._sample_aperture = sample_aperture_radius
        self._pixel_size_x = pixel_size_x
        self._pixel_size_y = pixel_size_y

        return

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


def calculate_momentum_transfer(ws):
    """Calculate momentum transfer including scaler momentum transfer Q, azimuthal projection (Qx) and Qy

    Note:
    - N: number of spectra
    - M: number of wave length bins
    - Qx: Azimuthal projection, as Q's projection to detector plane and then to X-axis (Q cos(phi))
    - Qy: Orthogonal to azimuthal projection, i.e., Q sin (phi)

    Parameters
    ----------
    ws : MatrixWorkspace
        Wavelength workspace to calculate momentum transfer from
    Returns
    -------
    ndarray, ndarray, ndarray, ndarray, ndarray
        N x M matrix for Q, N x M matrix for Qx, N x M matrix for Qy, (N, 1) array for 2theta,
        (N, 1) for sample pixel distance

    """
    # Check inputs and unit
    if ws is None:
        raise RuntimeError('Workspace cannot be None')
    elif ws.getAxis(0).getUnit().unitID() != 'Wavelength':
        raise RuntimeError('Input workspace {} for calculate Q resolution must be in unit Wavelength but not {}'
                           ''.format(ws, ws.getAxis(0).getUni().unitID()))

    # Obtain wave length
    wavelength_bin_boundary_matrix = ws.extractX()
    wavelength_bin_center_matrix = (wavelength_bin_boundary_matrix[:, 1:] +
                                    wavelength_bin_boundary_matrix[:, :-1]) / 2.0

    # Get instrument information
    spec_info = ws.spectrumInfo()
    two_theta_array = np.zeros((spec_info.size(), 1))  # 2-theta array is 2D (N, ) left to mask
    phi_array = np.zeros((spec_info.size(), 1))  # phi array is 2D (N x 1)
    sample_pixel_distance_array = np.zeros((spec_info.size(), 1))  # s2p array is 2D (N x 1)

    # Get instrument geometry information
    for i in range(spec_info.size()):
        if spec_info.hasDetectors(i) and not spec_info.isMonitor(i):
            # spectrum corresponds to a valid detector
            sample_pixel_distance_array[i] = spec_info.l2(i)
            two_theta_array[i] = spec_info.twoTheta(i)
            phi_array[i] = spec_info.azimuthal(i)
        else:
            # otherwise, we don't care
            two_theta_array[i] = np.nan
    # END-FOR

    # Calculate momentum transfer
    mask = np.isnan(two_theta_array.reshape((spec_info.size(),)))  # mask must use a (N,) array to act one a 2D
    q_array = 4.0 * np.pi * np.sin(0.5 * two_theta_array) / wavelength_bin_center_matrix
    q_array[mask] = 0.

    # Calculate Qx and Qy
    qx_array = np.cos(phi_array) * q_array
    qy_array = np.sin(phi_array) * q_array

    # Clean memory
    del mask, phi_array

    return q_array, qx_array, qy_array, two_theta_array, sample_pixel_distance_array


def dq2_geometry(L1, L2, R1, R2, wl, theta, pixel_size=0.007):
    r"""
    Geometry term for both dqx and dqy.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace
    L1: float
        source-to-sample distance (m)
    L2: float
        sample-to-detector distance (m)
    R1: float
        source aperture radius (m)
    R2: float
        sample aperture radius (m)
    wl: float
        wavelength mid-point (Angstrom)
    theta: float or ndarray(dtype=float)
        scattering angle (rad)
    pixel_size: float
        dimension of the pixel (m)

    Returns
    ------
    float
    """
    dq2 = 0.25 * np.square(L2 / L1 * R1)  \
        + 0.25 * np.square((L1 + L2) / L1 * R2) + np.square(pixel_size) / 12.0

    return dq2 * np.square(2.0 * np.pi * np.cos(theta)
                           * np.square(np.cos(2.0 * theta)) / wl / L2)


def dq2_gravity(L1, L2, wl, dwl, theta):
    r"""
    Gravity term for dqy.

    Parameters
    ----------
    L1: float
        source-to-sample distance (m)
    L2: float
        sample-to-detector distance (m)
    wl: float
        wavelength mid-point (Angstrom)
    dwl: float
        wavelength-spread (Angstrom)
    theta: float
        scattering angle (rad)

    Returns
    ------
    float
    """
    B = 0.5 * G_MN2_OVER_H2 * L2 * (L1 + L2)
    dq2 = 2. * np.square(B * wl * dwl) / 3.
    dq2 *= np.square(2.0 * np.pi * np.cos(theta)
                     * np.square(np.cos(2.0 * theta)) / wl / L2)
    # Converting from A^2 / m^4 to 1 / A^2
    return dq2 * 1.0e-40
