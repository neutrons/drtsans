# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
from drtsans.dataobjects import verify_same_q_bins
import numpy as np


__all__ = ['normalize_by_elastic_reference']


def normalize_by_elastic_reference(i_of_q, ref_i_of_q):

    # check i_of_q and ref_i_of_q shall have same binning
    if not verify_same_q_bins(i_of_q, ref_i_of_q):
        raise RuntimeError('Input I(Q) and elastic reference I(Q) have different Q and wavelength binning')

    # Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process
    q_min, q_max = determine_common_mod_q_range(i_of_q)

    # Find scale factor K(lambda) , that minimizes sum_q |refI(q, lambda_ref) - K(lambda) refI(q, lambda)|^2
    scale_factor_k_vec = calculate_scale_factor(ref_i_of_q, q_min, q_max)

    # Normalize input I(q, lambda)
    i_of_q = normalize_intensity(i_of_q, scale_factor_k_vec)

    return i_of_q, scale_factor_k_vec, sigma_scale_factor_k_vec


def determine_common_mod_q_range(iqmod):
    """Determine the common Q1D range among all the wavelengths such that I(q, lambda) does exist

    Detailed requirement:
        Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process

    Parameters
    ----------
    iqmod :  ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from

    Returns
    -------
    tuple
        minimum Q, maximum Q
    """
    # verify
    if iqmod.wavelength is None:
        raise RuntimeError('Input IQmod is required to have wavelength vector')

    # description of the algorithm for Q1D
    # get unique wave length
    wl_vec = np.unique(iqmod.wavelength)
    wl_vec.sort()

    # create a 2D array: wavelength, Q1D, intensity
    combo_vec = np.array([iqmod.wavelength, iqmod.mod_q, iqmod.intensity])
    # transpose for numpy selection and slicing
    combo_vec = combo_vec.transpose()

    # then for each unique wave length, determine the min and max
    min_q_vec = np.zeros_like(wl_vec, dtype=float)
    max_q_vec = np.zeros_like(wl_vec, dtype=float)
    for iw, wave_length in enumerate(wl_vec):
        # filter by wave length: column 0
        single_wl_q_vec = combo_vec[combo_vec[:, 0] == wave_length]
        # filter by intensity: column 2
        single_wl_q_vec = single_wl_q_vec[np.isfinite(single_wl_q_vec[:, 2])]
        # get minimum and max Q
        min_q_wl = np.min(single_wl_q_vec[:, 1])
        max_q_wl = np.max(single_wl_q_vec[:, 1])
        # set value
        min_q_vec[iw] = min_q_wl
        max_q_vec[iw] = max_q_wl

    # get the common qmin and qmax
    return np.max(min_q_vec), np.min(max_q_vec)


def calculate_scale_factor(i_of_q, q_min, q_max):
    """

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from
    q_min
    q_max

    Returns
    -------

    """
    # FIXME-Question: K is a function of wave length.  In principle, Q2D can result in same value?

    # Retrieve unique wave length in ascending order
    unique_wavelength_vec = np.unique(i_of_q.wavelength)
    assert len(unique_wavelength_vec.shape) == 0
    unique_wavelength_vec.sort()
    num_wavelength = unique_wavelength_vec.shape[0]

    # Retrieve reference wavelength for each q value within q_min and q_max
    ref_q_wl_vec = determine_reference_wavelength_q1d(i_of_q)

    # Create a matrix for q, wavelength, intensity and error
    # NOTE TODO : output I'(Q, lambda) will be a complete new i_of_q instance
    i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity,
                              i_of_q.error])
    i_q_wl_matrix.transpose()

    # Calculate P(wl), S(wl)
    p_vec = np.zeros_like(unique_wavelength_vec)
    s_vec = np.zeros_like(unique_wavelength_vec)
    k_vec = np.zeros_like(unique_wavelength_vec)
    k_error2_vec = np.zeros_like(unique_wavelength_vec)

    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # filter i_q_wl_matrix for certain wavelength
        i_q_matrix = i_q_wl_matrix[i_q_wl_matrix[:, 1] == lambda_i]
        # get the index of q_min and q_max
        # TODO - need to test
        i_q_min = np.argmin(i_q_matrix[:, 0], q_min)[0]
        i_q_max = np.argmax(i_q_matrix[:, 0], q_max)[1]
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            # FIXME -  A better data structure is required for this and the following is wrong!
            ref_wl = ref_q_wl_vec[q_i]
            p_vec[i_wl] += get_intensity(i_of_q, q_i, ref_wl) * i_q_matrix[i_q][2]
            s_vec[i_wl] += i_q_matrix[i_q][2] * i_q_matrix[i_q][2]

        # calculate K(wl)
        k_vec[i_wl] = p_vec[i_wl] / s_vec[i_wl]

        # calculate delta K(wl)
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            ref_wl = ref_q_wl_vec[q_i]
            # delta I(q, lambda) = delta I^{lambda}(q)
            term0 = i_q_matrix[i_q][3]
            # ... ...
            term1 = (get_intensity(i_of_q, q_i, ref_wl)**2 - 2 * i_q_matrix[q_i][2] * p_vec[i_wl])/s_vec[i_wl]**2
            # delta I(q, lambda^ref)
            term2 = get_error(i_of_q, q_i, ref_wl)
            # I(q, lambda)/S(lambda) = I^{lambda}(q) / S(lambda)
            term3 = i_q_matrix[q_i][2] / s_vec[q_i]
            k_error2_vec[i_wl] += (term0 * term1) ** 2 + (term2 * term3)**2

    return k_vec, k_error2_vec, p_vec, s_vec


def calculate_scale_factor_mesh_grid(i_of_q, q_min, q_max):
    """

    Same functionality as calculate_scale_factor but the algorithm is improved
    as I(Q, mavelength) are in meshgrid

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from
    q_min
    q_max

    Returns
    -------

    """
    return [], [], [], []


def determine_reference_wavelength_q1d(i_of_q):
    """Determine the reference wavelength for each Q.

    The reference wavelength of a specific Q or (qx, qy)
    is defined as the shortest wavelength for all the finite I(Q, wavelength) or
    I(qx, qy, wavelength)

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        I(Q, wavelength)

    Returns
    -------
    numpy.ndarray
        2D array of (Q, wavelength)

    """
    # Construct nd array for Q, I and wavelength
    combo_matrix = np.array([i_of_q.mod_q, i_of_q.intensity, i_of_q.wavelength])
    combo_matrix = combo_matrix.transpose()  # good for filter and slicing
    # Create a list of unique Q
    unique_q_vec = np.unique(i_of_q.mod_q)

    # Init return vector
    ref_wavelength_vec = np.ndarray(shape=(unique_q_vec.shape[0], 2), dtype='float')

    # Remove all the I(q, wl) with intensity as nan or infinity:

    # For each Q, search wavelength min
    for index, q_value in enumerate(unique_q_vec):
        # filter out the items with desired Q value
        filtered_combo_matrix = combo_matrix[combo_matrix[:, 0] == q_value]
        filtered_combo_matrix = filtered_combo_matrix[np.isfinite(filtered_combo_matrix[:, 1])]
        # find minimum wavelength and add to output array
        min_wl = np.min(filtered_combo_matrix[:, 2])
        # set
        ref_wavelength_vec[index][0] = q_value
        ref_wavelength_vec[index][1] = min_wl

    return ref_wavelength_vec


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


def normalize_intensity(i_of_q, k_vec, k_error_vec, ref_wl_vec, p_vec, s_vec,
                        unique_wavelength_vec):
    """

    Parameters
    ----------
    i_of_q
    k_vec
    k_error_vec
    ref_wl_vec
    p_vec
    s_vec

    Returns
    -------

    """
    # Create a matrix for q, wavelength, intensity and error
    # NOTE TODO : output I'(Q, lambda) will be a complete new i_of_q instance
    i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity,
                              i_of_q.error])
    i_q_wl_matrix.transpose()

    # Output

    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # filter i_q_wl_matrix for certain wavelength
        i_q_matrix = i_q_wl_matrix[i_q_wl_matrix[:, 1] == lambda_i]

        # normalize intensity by K^{lambda} * I^{lambda}(q)
        corrected_intensity_vec = i_q_matrix[:, 2] * k_vec[i_wl]

        # new matrix
        new_matrix = np.array([corrected_intensity_vec])

    # construct a new instance of I(Q)
    full_matrix = concatenate_matrix(new_matrix)
    corrected_i_of_q = IQmod(full_matri[0], full_matrix[1], i_of_q.mod_q, i_of_q.delta_q, i_of_q.wavelength)
    
    return i_of_q








