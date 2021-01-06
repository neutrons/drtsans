# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
from drtsans.dataobjects import verify_same_q_bins
import numpy as np
from dataclasses import dataclass


__all__ = ['normalize_by_elastic_reference']


@dataclass
class ReferenceWavelengths:
    """
    Class for keeping track of reference wavelength for each momentum transfer Q (1D)
    """
    q_vec: np.ndarray
    ref_wl_vec: np.ndarray
    intensity_vec: np.ndarray
    error_vec: np.ndarray

    def __init__(self, q_values, ref_wavelengths, intensities, errors):
        """

        Parameters
        ----------
        q_values
        ref_wavelengths
        intensities
        errors
        """
        self.q_vec = q_values
        self.ref_wl_vec = ref_wavelengths
        self.intensity_vec = intensities
        self.error_vec = errors


def reshape_q_wavelength_matrix(i_of_q):
    """Reshape I(Q) into a mesh grid of (Q, wavelength) and limit Q into q_min and q_max

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from

    Returns
    -------
    tuple
        wavelength vector, q vector,  intensity (2D), error (2D)

    """
    # Retrieve unique wave length in ascending order
    wavelength_vec = np.unique(i_of_q.wavelength)
    assert len(wavelength_vec.shape) == 1
    wavelength_vec.sort()

    q_vec = np.unique(i_of_q.mod_q)
    assert len(q_vec.shape) == 1
    q_vec.sort()

    # Create a matrix for q, wavelength, intensity and error
    i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity,
                              i_of_q.error])
    i_q_wl_matrix = i_q_wl_matrix.transpose()

    # Order by wavelength and momentum transfer (Q)
    i_q_wl_matrix = i_q_wl_matrix[np.lexsort((i_q_wl_matrix[:, 1], i_q_wl_matrix[:, 0]))]

    # Unique wavelength and unique momentum transfer
    wl_vector = np.unique(i_of_q.wavelength)
    q_vector = np.unique(i_of_q.mod_q)
    # verify whether (q, wl) are on mesh grid by checking unique Q and wavelength
    assert wl_vector.shape[0] * q_vector.shape[0] == i_of_q.intensity.shape[0]

    # Reformat
    intensity_array = i_q_wl_matrix[:, 2].reshape((q_vector.shape[0], wl_vector.shape[0]))
    error_array = i_q_wl_matrix[:, 3].reshape((q_vector.shape[0], wl_vector.shape[0]))

    return wl_vector, q_vector, intensity_array, error_array


def normalize_by_elastic_reference(i_of_q, ref_i_of_q):

    # check i_of_q and ref_i_of_q shall have same binning
    if not verify_same_q_bins(i_of_q, ref_i_of_q):
        raise RuntimeError('Input I(Q) and elastic reference I(Q) have different Q and wavelength binning')

    # Reshape Q, wavelength, intensities and errors to unique 1D array or 2D array
    wl_vec, q_vec, i_array, error_array = reshape_q_wavelength_matrix(i_of_q)

    # Calculate Qmin and Qmax
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # Calculate reference
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array,
                                                        qmin_index, qmax_index)

    # Calculate scale factor
    k_vec, k_error_vec, p_vec, s_vec = calculate_scale_factor_mesh_grid(wl_vec, i_array, error_array,
                                                                        ref_wl_ie, qmin_index, qmax_index)

    # Normalize
    normalized = normalize_intensity_q1d(wl_vec, q_vec, i_array, error_array,
                                         ref_wl_ie, k_vec, p_vec, s_vec,
                                         qmin_index, qmax_index)

    # TODO - convert normalized intensities and errors to ImodQ

    return normalized, k_vec, k_error_vec


def determine_common_mod_q_range_mesh(q_vec, intensity_array):
    """Determine the common Q1D range among all the wavelengths such that I(q, lambda) does exist.

    This method assumes that I(Q, wavelength) are on mesh grid of Q and wavelength

    Detailed requirement:
        Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process

    Parameters
    ----------
    q_vec: numpy.ndarray
        vector of sorted unique Q
    intensity_array: numpy.ndarray
        2D array of intensity.  Each row is of same wavelength

    Returns
    -------
    tuple
        index of qmin and qmax

    """
    # Find q min
    qmin_index = None
    qmax_index = None

    # Sanity check
    assert q_vec.shape[0] == intensity_array.shape[0], 'Shape mismatch'

    num_q = q_vec.shape[0]
    for q_index in range(num_q):
        if len(np.where(np.isnan(intensity_array[q_index]))[0]) == 0:
            qmin_index = q_index
            break
    for q_index in range(num_q - 1, -1, -1):
        if len(np.where(np.isnan(intensity_array[q_index]))[0]) == 0:
            qmax_index = q_index
            break

    if qmin_index is None:
        raise RuntimeError('Unable to find common q range')

    return qmin_index, qmax_index


def calculate_scale_factor_mesh_grid(wl_vec, intensity_array, error_array,
                                     ref_wl_intensities, qmin_index, qmax_index):
    """Same functionality as calculate_scale_factor but the algorithm is improved
    as I(Q, wavelength) are in meshgrid

    Parameters
    ----------
    wl_vec: numpy.array
        wavelength vector
    intensity_array: numpy.array
        intensity 2D array
    error_array: numpy.array
        error 2D array
    ref_wl_intensities: ReferenceWavelengths
        reference wavelength intensity/error
    qmin_index: int
        index of min Q in q vector
    qmax_index: int
        index of max Q in q vector

    Returns
    -------
    tuple
        K vector, K error vector, P vector, S vector
    """
    # Check input
    assert wl_vec.shape[0] == intensity_array.shape[1]

    # Calculate P(wl), S(wl)
    p_vec = np.zeros_like(wl_vec)
    s_vec = np.zeros_like(wl_vec)
    k_error2_vec = np.zeros_like(wl_vec)

    for i_wl, lambda_i in enumerate(wl_vec):
        # P(wl) = sum_q I(q, ref_wl) * I(q, wl)
        p_value = np.sum(ref_wl_intensities.intensity_vec[qmin_index:qmax_index + 1] *
                         intensity_array[:, i_wl][qmin_index:qmax_index + 1])
        # S(wl) = sum_q I(q, wl)**2
        s_value = np.sum(intensity_array[:, i_wl][qmin_index:qmax_index + 1]**2)

        # assign
        p_vec[i_wl] = p_value
        s_vec[i_wl] = s_value

        term0 = error_array[:, i_wl][qmin_index:qmax_index+1]
        term1 = (ref_wl_intensities.intensity_vec[qmin_index:qmax_index+1] * s_value -
                 2. * intensity_array[:, i_wl][qmin_index:qmax_index+1] * p_value) / s_value**2
        term2 = ref_wl_intensities.error_vec[qmin_index:qmax_index+1]
        term3 = intensity_array[:, i_wl][qmin_index:qmax_index+1] / s_value

        k_error2_vec[i_wl] = np.sum((term0 * term1)**2 + (term2 * term3)**2)

    # Calculate K
    k_vec = p_vec / s_vec

    return k_vec, np.sqrt(k_error2_vec), p_vec, s_vec


def determine_reference_wavelength_q1d_mesh(wavelength_vec, q_vec, intensity_array, error_array,
                                            qmin_index, qmax_index):
    """Determine the reference wavelength for each Q.

    The reference wavelength of a specific Q or (qx, qy)
    is defined as the shortest wavelength for all the finite I(Q, wavelength) or
    I(qx, qy, wavelength)

    Parameters
    ----------
    wavelength_vec: numpy.ndarray
        ...
    q_vec: numpy.ndarray
        ...
    intensity_array: numpy.ndarray
        ...
    error_array: numpy.ndarray
        ...
    qmin_index: int
        index of qmin in q-vector
    qmax_index: int
        index of qmax in q-vector

    Returns
    -------
    ReferenceWavelengths
        Reference wavelengths for each momentum transfer Q and the corresponding intensity and error

    """
    # Sanity check
    assert wavelength_vec.shape[0] == intensity_array.shape[1]

    # Minimum wavelength bin is the reference wavelength
    min_wl_vec = np.zeros_like(q_vec) + wavelength_vec[0]

    # Minimum intensity and error
    min_intensity_vec = np.copy(intensity_array[:, 0])
    min_error_vec = np.copy(error_array[:, 0])

    # Set the unused defined reference wavelength (outside of qmin and qmax)'s
    # intensity and error to nan
    min_intensity_vec[0:qmin_index] = np.nan
    min_intensity_vec[qmax_index+1:] = np.nan
    min_error_vec[0:qmin_index] = np.nan
    min_error_vec[qmax_index+1:] = np.nan

    return ReferenceWavelengths(q_vec, min_wl_vec, min_intensity_vec, min_error_vec)


def determine_reference_wavelength_q2d(i_of_q):
    """Determine the reference wavelength for each Q.

    The reference wavelength of a specific Q or (qx, qy)
    is defined as the shortest wavelength for all the finite I(Q, wavelength) or
    I(qx, qy, wavelength)

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        I(qx, qy, wavelength)

    Returns
    -------
    numpy.ndarray
        3D array of (qx, qy, wavelength)

    """
    # Construct nd array for qx, qy, I and wavelength
    combo_matrix = np.array([i_of_q.qx, i_of_q.qy, i_of_q.intensity, i_of_q.wavelength])
    combo_matrix = combo_matrix.transpose()  # good for filter and slicing
    # Create a list of unique Q
    q2d_vec = np.array([i_of_q.qx, i_of_q.qy])
    q2d_vec.transpose()
    unique_q_vec = np.unique(q2d_vec, axis=0)

    # Init return vector
    ref_wavelength_vec = np.ndarray(shape=(unique_q_vec.shape[0], 3), dtype='float')

    # Remove all the I(q, wl) with intensity as nan or infinity
    combo_matrix = combo_matrix[np.isfinite(combo_matrix)]

    # For each Q, search wavelength min
    for index, q_value in unique_q_vec:
        # filter out the items with desired Q value
        combo_filter_qx = combo_matrix[combo_matrix[:, 0] == q_value[0]]
        combo_filter_qy = combo_filter_qx[combo_filter_qx[:, 1] == q_value[1]]
        # find minimum wavelength and add to output array
        min_wl = np.min(combo_filter_qy[:, 2])
        # set value
        ref_wavelength_vec[index][0] = q_value[0]
        ref_wavelength_vec[index][1] = q_value[1]
        ref_wavelength_vec[index][2] = min_wl

    return ref_wavelength_vec


def normalize_intensity_q1d(wl_vec, q_vec, intensity_array, error_array, ref_wl_ints_errs, k_vec, p_vec, s_vec,
                            qmin_index, qmax_index):
    """Normalize Q1D intensities and errors

    Parameters
    ----------
    wl_vec
    q_vec
    intensity_array
    error_array
    ref_wl_ints_errs
    k_vec
    p_vec
    s_vec
    qmin_index
    qmax_index

    Returns
    -------

    """

    # Sanity check
    assert wl_vec.shape[0] == intensity_array.shape[1]
    assert q_vec.shape[0] == error_array.shape[0]
    assert intensity_array.shape == error_array.shape

    # Normalized intensities
    normalized_intensity_array = intensity_array * k_vec
    normalized_error2_array = np.zeros_like(error_array)

    # Lowest wavelength bin does not require normalization as K = 1, i_wl = 0
    normalized_error2_array[:, 0] = error_array[:, 0]**2

    # Reshape
    ri_vec = ref_wl_ints_errs.intensity_vec.reshape((q_vec.shape[0], 1))
    re_vec = ref_wl_ints_errs.error_vec

    # qmax is included.  need i_qmax to slicing
    i_qmax = qmax_index + 1

    # Loop over wavelength
    num_wl = wl_vec.shape[0]
    for i_wl in range(1, num_wl):

        intensity_vec = intensity_array[:, i_wl].reshape((q_vec.shape[0], 1))

        # Calculate Y: Y_ij = I_i * R_j * s - I_i * 2 * I_j * p
        y_matrix = \
            intensity_vec * (ri_vec.transpose()) * s_vec[i_wl] - \
            intensity_vec * (intensity_vec.transpose()) * (2 * p_vec[i_wl])
        y_diag = np.diag(y_matrix)
        # y_matrix[i, :] corresponds to a single q_i/r_i
        # y_matrix[:, j] corresponds to a single q_j/r_j

        # Term 2
        # t2 += [delta I(q', wl)]**2 * Y(q, q'', wl)**2 / S(lw)**4
        t2sum_vec = \
            error_array[qmin_index:i_qmax, i_wl]**2 * y_matrix[:, qmin_index:i_qmax]**2 / s_vec[i_wl]**4

        # Term 3
        # t3 += [delta I(q_j, ref_wl[q_j]]^2 * [I(q_j, wl) * I(q, wl)]^2 / S(wl)^2
        t3sum_vec = \
            intensity_array[:, i_wl]**2 * np.sum(
                re_vec[qmin_index:i_qmax]**2 * intensity_array[qmin_index:i_qmax, i_wl]**2 / s_vec[i_wl]**2)

        # term 1
        # outside of qmin and qmax: t1 = [delta I(q, wl)]**2 * [P(wl) / S(wl)]**2
        t1sum_vec = (error_array[:, i_wl] * p_vec[i_wl] / s_vec[i_wl])**2
        t1sum_vec[qmin_index:qmax_index+1] =\
            error_array[qmin_index:qmax_index+1, i_wl]**2\
            * (p_vec[i_wl]**2 * s_vec[i_wl]**2
               + 2 * p_vec[i_wl] * s_vec[i_wl] * y_diag[qmin_index:qmax_index+1]) / s_vec[i_wl]**4

        # sum up
        normalized_error2_array[:, i_wl] += t1sum_vec + t2sum_vec.sum(axis=1) + t3sum_vec

    return normalized_intensity_array, np.sqrt(normalized_error2_array)
