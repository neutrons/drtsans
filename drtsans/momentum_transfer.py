"""
Resolution calculations common to all SANS
"""
import numpy as np
from numpy import linalg
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
    """Calculate momentum transfer in vector by operating on source, sample and detector positions,
    i.e., vec{q} = 2 pi (vec{k_out} - vec{k_in}) / lambda

    Note:
    - N: number of spectra
    - M: number of wave length bins
    - Qx: Q component at X-Z plane
    - Qy: Q component at Y plane

    Parameters
    ----------
    ws : MatrixWorkspace
        Wavelength workspace to calculate momentum transfer from
    Returns
    -------
    ndarray, ndarray, ndarray, ndarray
        N x M matrix for Q, N x M matrix for Qx, N x M matrix for Qy, N array for 2theta,

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
    num_spec = spec_info.size()
    two_theta_array = np.zeros(spec_info.size())

    # sample and moderator information: get K_i
    sample_pos = ws.getInstrument().getSample().getPos()
    source_pos = ws.getInstrument().getSource().getPos()
    k_in = sample_pos - source_pos
    k_in /= linalg.norm(k_in)

    unit_q_vector = np.zeros(shape=(num_spec, 3), dtype='float')

    # TODO - k_out can be calculated in the vector operation by recording det_i_pos, the l2 is almost "free"
    for iws in range(num_spec):
        if spec_info.hasDetectors(iws) and not spec_info.isMonitor(iws):
            # calculate Qx, Qy, Qz
            det_i_pos = ws.getDetector(iws).getPos()
            k_out = det_i_pos - sample_pos
            k_out /= linalg.norm(k_out)
            unit_q_vector[iws] = k_out - k_in

            # 2theta
            two_theta_array[iws] = spec_info.twoTheta(iws)
        # otherwise, unit_q_vector[iws] is zero
    # END-FOR

    # Calculate Q from unit Q vector and 2pi/lambda
    qx_matrix = 2.*np.pi*unit_q_vector[:, 0].reshape((num_spec, 1)) / wavelength_bin_center_matrix
    qy_matrix = 2.*np.pi*unit_q_vector[:, 1].reshape((num_spec, 1)) / wavelength_bin_center_matrix
    qz_matrix = 2.*np.pi*unit_q_vector[:, 2].reshape((num_spec, 1)) / wavelength_bin_center_matrix

    q_matrix = np.sqrt(qx_matrix**2 + qy_matrix**2 + qz_matrix**2)
    q_theta_matrix = np.sqrt(qx_matrix**2 + qz_matrix**2)

    # print('[DEBUG] Q  matrix: shape={}\n{}'.format(q_matrix.shape, q_matrix))
    # print('[DEBUG] Qx matrix: shape={}\n{}'.format(qx_matrix.shape, qx_matrix))
    # print('[DEBUG] Qy matrix: shape={}\n{}'.format(qy_matrix.shape, qy_matrix))
    # print('[DEBUG] Qz matrix: shape={}\n{}'.format(qz_matrix.shape, qz_matrix))

    return q_matrix, q_theta_matrix, qy_matrix, two_theta_array


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
    theta: float
        scattering angle (rad)
    pixel_size: float
        dimension of the pixel (m)

    Returns
    ------
    float
    """
    dq2 = 0.25 * np.square(L2 / L1 * R1)  \
        + 0.25 * np.square((L1 + L2) / L1 * R2) + np.square(pixel_size) / 12.0

    print('[DEBUG RES] factor1 = {},  dq2 = {}'
          ''.format(np.square(2.0 * np.pi * np.cos(theta)
                           * np.square(np.cos(2.0 * theta)) / wl / L2),  dq2))

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
    print('[RICARDO] B = {}, Gravity = {}'.format(B, dq2))
    dq2 *= np.square(2.0 * np.pi * np.cos(theta)
                     * np.square(np.cos(2.0 * theta)) / wl / L2)
    # Converting from A^2 / m^4 to 1 / A^2
    return dq2 * 1.0e-40
