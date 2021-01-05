# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
from drtsans.dataobjects import verify_same_q_bins
from drtsans.dataobjects import IQmod
import numpy as np


__all__ = ['normalize_by_elastic_reference']


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

    # Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process
    q_min, q_max = determine_common_mod_q_range(i_of_q)

    # Find scale factor K(lambda) , that minimizes sum_q |refI(q, lambda_ref) - K(lambda) refI(q, lambda)|^2
    array_tuples = calculate_scale_factor(ref_i_of_q, q_min, q_max)
    scale_factor_k_vec = array_tuples[0]
    delta_k_vec = array_tuples[1]
    p_factor_vec = array_tuples[2]
    s_factor_vec = array_tuples[3]
    unique_wl_vec = array_tuples[4]
    ref_wl_vec = array_tuples[5]

    # Normalize input I(q, lambda)
    i_of_q = normalize_intensity(i_of_q, scale_factor_k_vec, ref_wl_vec,
                                 p_factor_vec, s_factor_vec, unique_wl_vec, q_min, q_max)

    return i_of_q, scale_factor_k_vec, delta_k_vec


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
    tuple
        K, delta K, P, S, unique Q, reference wavelength I(q) and dI(q)

    """
    # FIXME-Question: K is a function of wave length.  In principle, Q2D can result in same value?

    # Retrieve unique wave length in ascending order
    unique_wavelength_vec = np.unique(i_of_q.wavelength)
    assert len(unique_wavelength_vec.shape) == 1
    unique_wavelength_vec.sort()

    # Retrieve reference wavelength for each q value within q_min and q_max
    ref_q_wl_vec = determine_reference_wavelength_q1d(i_of_q)

    # Create a matrix for q, wavelength, intensity and error
    # NOTE TODO : output I'(Q, lambda) will be a complete new i_of_q instance
    i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity,
                              i_of_q.error])
    i_q_wl_matrix = i_q_wl_matrix.transpose()
    # print(f'Q-WL-I-Sigma matrix shape = {i_q_wl_matrix.shape}')

    # Calculate P(wl), S(wl)
    p_vec = np.zeros_like(unique_wavelength_vec)
    s_vec = np.zeros_like(unique_wavelength_vec)
    k_vec = np.zeros_like(unique_wavelength_vec)
    k_error2_vec = np.zeros_like(unique_wavelength_vec)

    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # filter i_q_wl_matrix for certain wavelength
        i_q_matrix = i_q_wl_matrix[i_q_wl_matrix[:, 1] == lambda_i]
        # print(f'[DEBUG-INFO] wavelength = {lambda_i}, values size = {i_q_matrix.shape}')

        # get the index of q_min and q_max
        i_q_min = np.argmin(np.abs(i_q_matrix[:, 0] - q_min))  # numpy int64
        assert q_min == i_q_matrix[:, 0][i_q_min]
        i_q_max = np.argmin(np.abs(i_q_matrix[:, 0] - q_max))
        assert q_max == i_q_matrix[:, 0][i_q_max]

        # loop over each Q for S(wl) and P(wl)
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            # print(f'Over {i_q}-th Q {q_i}')
            # acquire index to reference lambda
            ref_index = np.argmin(np.abs(ref_q_wl_vec[:, 0] - q_i))
            # ref_wl = ref_q_wl_vec[ref_index, 1]
            i_q_ref_wl = ref_q_wl_vec[ref_index, 2]
            # print(f'  ref index = {ref_index}  '
            #       f'reference wl = {ref_wl}  I(q, ref) = {i_q_ref_wl}  I(q, wl) = {i_q_matrix[i_q][2]}')
            p_vec[i_wl] += i_q_ref_wl * i_q_matrix[i_q][2]
            # print(f'q-index {i_q}  increment P = {i_q_ref_wl * i_q_matrix[i_q][2]}')
            s_vec[i_wl] += i_q_matrix[i_q][2] * i_q_matrix[i_q][2]

        # calculate K(wl)
        k_vec[i_wl] = p_vec[i_wl] / s_vec[i_wl]

        # calculate delta K(wl)
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            # acquire index to reference lambda
            ref_index = np.argmin(np.abs(ref_q_wl_vec[:, 0] - q_i))
            i_q_ref_wl = ref_q_wl_vec[ref_index, 2]
            err_q_ref_wl = ref_q_wl_vec[ref_index, 3]
            # delta I(q, lambda) = delta I^{lambda}(q)
            term0 = i_q_matrix[i_q, 3]
            # I(q, ref_wl(q)) * S(wl) - 2 * I(q, wl) * P(wl) / S(wl)**2
            term1 = (i_q_ref_wl * s_vec[i_wl] - 2 * i_q_matrix[i_q][2] * p_vec[i_wl])/s_vec[i_wl]**2
            # delta I(q, lambda^ref)
            term2 = err_q_ref_wl
            # I(q, lambda)/S(lambda) = I^{lambda}(q) / S(lambda)
            term3 = i_q_matrix[i_q, 2] / s_vec[i_wl]
            # increment = (t0 * t1)**2 + (t2 * t3)**2
            k_error2_vec[i_wl] += (term0 * term1) ** 2 + (term2 * term3)**2

            # print(f'  error(q, wl)    t0 = {term0}')
            # print(f'                  t1 = {term1}')
            # print(f'                        I(Q, refWL) = {i_q_ref_wl}')
            # print(f'                        S(wl      ) = {s_vec[i_wl]}')
            # print(f'                        I(Q,    wl) = {i_q_matrix[i_q][2]}')
            # print(f'                        P(wl      ) = {p_vec[i_wl]}')
            # print(f'  reference error t2 = {term2}')
            # print(f'  t3 = {term3}')
            # print(f'  Increment = {(term0 * term1) ** 2 + (term2 * term3)**2}')

    # END-FOR

    # Get K error vector
    k_error_vec = np.sqrt(k_error2_vec)

    return k_vec, k_error_vec, p_vec, s_vec, unique_wavelength_vec, ref_q_wl_vec


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
    # NOTE: output I'(Q, lambda) will be a complete new i_of_q instance


    # Limit the wavelength between qmin and qmax
    qmin_index = np.argmin(np.abs(unique_wavelength_vec - q_min))
    qmax_index = np.argmin(np.abs(unique_wavelength_vec - q_max))
    print(f'Qmin @ {qmin_index}, Qmax @ {qmax_index}')

    reshape_q_wavelength_matrix(i_of_q, q_min, q_max)
    raise RuntimeError('... ...')

    # print(f'Q-WL-I-Sigma matrix shape = {i_q_wl_matrix.shape}')

    # Calculate P(wl), S(wl)
    p_vec = np.zeros_like(unique_wavelength_vec)
    s_vec = np.zeros_like(unique_wavelength_vec)
    k_vec = np.zeros_like(unique_wavelength_vec)
    k_error2_vec = np.zeros_like(unique_wavelength_vec)

    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # filter i_q_wl_matrix for certain wavelength
        i_q_matrix = i_q_wl_matrix[i_q_wl_matrix[:, 1] == lambda_i]
        # print(f'[DEBUG-INFO] wavelength = {lambda_i}, values size = {i_q_matrix.shape}')

        # get the index of q_min and q_max
        i_q_min = np.argmin(np.abs(i_q_matrix[:, 0] - q_min))  # numpy int64
        assert q_min == i_q_matrix[:, 0][i_q_min]
        i_q_max = np.argmin(np.abs(i_q_matrix[:, 0] - q_max))
        assert q_max == i_q_matrix[:, 0][i_q_max]

        # loop over each Q for S(wl) and P(wl)
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            # print(f'Over {i_q}-th Q {q_i}')
            # acquire index to reference lambda
            ref_index = np.argmin(np.abs(ref_q_wl_vec[:, 0] - q_i))
            # ref_wl = ref_q_wl_vec[ref_index, 1]
            i_q_ref_wl = ref_q_wl_vec[ref_index, 2]
            # print(f'  ref index = {ref_index}  '
            #       f'reference wl = {ref_wl}  I(q, ref) = {i_q_ref_wl}  I(q, wl) = {i_q_matrix[i_q][2]}')
            p_vec[i_wl] += i_q_ref_wl * i_q_matrix[i_q][2]
            # print(f'q-index {i_q}  increment P = {i_q_ref_wl * i_q_matrix[i_q][2]}')
            s_vec[i_wl] += i_q_matrix[i_q][2] * i_q_matrix[i_q][2]

        # calculate K(wl)
        k_vec[i_wl] = p_vec[i_wl] / s_vec[i_wl]

        # calculate delta K(wl)
        for i_q in range(i_q_min, i_q_max+1):
            q_i = i_q_matrix[:, 0][i_q]
            # acquire index to reference lambda
            ref_index = np.argmin(np.abs(ref_q_wl_vec[:, 0] - q_i))
            i_q_ref_wl = ref_q_wl_vec[ref_index, 2]
            err_q_ref_wl = ref_q_wl_vec[ref_index, 3]
            # delta I(q, lambda) = delta I^{lambda}(q)
            term0 = i_q_matrix[i_q, 3]
            # I(q, ref_wl(q)) * S(wl) - 2 * I(q, wl) * P(wl) / S(wl)**2
            term1 = (i_q_ref_wl * s_vec[i_wl] - 2 * i_q_matrix[i_q][2] * p_vec[i_wl])/s_vec[i_wl]**2
            # delta I(q, lambda^ref)
            term2 = err_q_ref_wl
            # I(q, lambda)/S(lambda) = I^{lambda}(q) / S(lambda)
            term3 = i_q_matrix[i_q, 2] / s_vec[i_wl]
            # increment = (t0 * t1)**2 + (t2 * t3)**2
            k_error2_vec[i_wl] += (term0 * term1) ** 2 + (term2 * term3)**2

            # print(f'  error(q, wl)    t0 = {term0}')
            # print(f'                  t1 = {term1}')
            # print(f'                        I(Q, refWL) = {i_q_ref_wl}')
            # print(f'                        S(wl      ) = {s_vec[i_wl]}')
            # print(f'                        I(Q,    wl) = {i_q_matrix[i_q][2]}')
            # print(f'                        P(wl      ) = {p_vec[i_wl]}')
            # print(f'  reference error t2 = {term2}')
            # print(f'  t3 = {term3}')
            # print(f'  Increment = {(term0 * term1) ** 2 + (term2 * term3)**2}')

    # END-FOR

    # Get K error vector
    k_error_vec = np.sqrt(k_error2_vec)

    return k_vec, k_error_vec, p_vec, s_vec, unique_wavelength_vec, ref_q_wl_vec


def determine_reference_wavelength_q1d(i_of_q):
    """Determine the reference wavelength for each Q.

    The reference wavelength of a specific Q or (qx, qy)
    is defined as the shortest wavelength for all the finite I(Q, wavelength) or
    I(qx, qy, wavelength)

    They shall be the same between qmin and qmax

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        I(Q, wavelength)

    Returns
    -------
    numpy.ndarray
        2D array of (Q, wavelength, intensity, error)

    """
    # Construct nd array for Q, I and wavelength
    combo_matrix = np.array([i_of_q.mod_q, i_of_q.intensity, i_of_q.wavelength, i_of_q.error])
    combo_matrix = combo_matrix.transpose()  # good for filter and slicing
    # Create a list of unique Q
    unique_q_vec = np.unique(i_of_q.mod_q)

    # Init return vector
    ref_wavelength_vec = np.ndarray(shape=(unique_q_vec.shape[0], 4), dtype='float')

    # Remove all the I(q, wl) with intensity as nan or infinity:

    # For each Q, search wavelength min
    for index, q_value in enumerate(unique_q_vec):
        # filter out the items with desired Q value
        filtered_combo_matrix = combo_matrix[combo_matrix[:, 0] == q_value]
        filtered_combo_matrix = filtered_combo_matrix[np.isfinite(filtered_combo_matrix[:, 1])]
        # find minimum wavelength and add to output array
        min_index = np.argmin(filtered_combo_matrix[:, 2])
        # set values
        ref_wavelength_vec[index][0] = q_value
        ref_wavelength_vec[index][1] = filtered_combo_matrix[:, 2][min_index]
        ref_wavelength_vec[index][2] = filtered_combo_matrix[:, 1][min_index]
        ref_wavelength_vec[index][3] = filtered_combo_matrix[:, 3][min_index]

    return ref_wavelength_vec


from dataclasses import dataclass

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


def determine_reference_wavelength_q1d_mesh(wavelength_vec, q_vec, intensity_array, error_array):
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

    Returns
    -------
    ReferenceWavelengths
        Reference wavelengths for each momentum transfer Q and the corresponding intensity and error

    """
    # Sanity check
    assert wavelength_vec.shape[0] == intensity_array.shape[1]

    # Get the maximum value of intensity
    max_intensity = intensity_array[np.isfinite(intensity_array)].max()

    # Replace nan by super max value
    super_max_val = max_intensity + 3.
    intensity_array[np.isnan(intensity_array)] = super_max_val

    # Get the minimum index of each wavelength
    min_index_vec = np.argmin(intensity_array, axis=1)
    min_wl_vec = wavelength_vec[min_index_vec]

    # Recover nan values for intensities
    intensity_array[intensity_array > max_intensity + 1.] = np.nan

    # Retrieve minimum wavelength for each

    # Get the position arrays for 2D matrix and then slice the intensity and error
    where_tuple = (np.arange(q_vec.shape[0]), min_index_vec)

    min_intensity_vec = intensity_array[where_tuple]
    min_error_vec = intensity_array[where_tuple]

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


def normalize_intensity(i_of_q, k_vec, ref_wl_vec, p_vec, s_vec,
                        unique_wavelength_vec, q_min, q_max):
    """

    Parameters
    ----------
    i_of_q
    k_vec
    ref_wl_vec
    p_vec
    s_vec
    unique_wavelength_vec

    Returns
    -------

    """
    # Create a matrix for q, wavelength, intensity and error
    # output I'(Q, lambda) will be a complete new i_of_q instance
    if i_of_q.delta_mod_q is None:
        # create a fake deltaQ for slicing
        delta_q_vec = i_of_q.mod_q
        delta_q_none = True
    else:
        delta_q_vec = i_of_q.delta_mod_q
        delta_q_none = False

    # TODO FIXME - this can be improved/refactor under the assumption that (Q, wavelength) are on mesh grid
    # unique Q vector
    unique_q_vec = np.unique(i_of_q.mod_q)
    unique_q_vec.sort()
    qmin_index = np.argmin(np.abs(unique_q_vec - q_min))
    qmax_index = np.argmin(np.abs(unique_q_vec - q_max))
    print(f'qmin = {q_min} @ {qmin_index}; qmax = {q_max} @ {qmax_index}')

    i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity,
                              i_of_q.error, delta_q_vec])
    i_q_wl_matrix = i_q_wl_matrix.transpose()

    # Init arrays to start
    new_mod_q = np.ndarray(shape=(0,), dtype=i_q_wl_matrix.dtype)
    new_q_error = np.ndarray(shape=(0,), dtype=i_q_wl_matrix.dtype)
    new_wavelength = np.ndarray(shape=(0,), dtype=i_q_wl_matrix.dtype)
    new_intensity = np.ndarray(shape=(0,), dtype=i_q_wl_matrix.dtype)
    new_error_sq = np.ndarray(shape=(0,), dtype=i_q_wl_matrix.dtype)

    # Calculate the normalized intensity error
    for i_wl, lambda_i in enumerate(unique_wavelength_vec):
        # filter i_q_wl_matrix for certain wavelength
        i_q_matrix = i_q_wl_matrix[i_q_wl_matrix[:, 1] == lambda_i]
        print(f'[DEBUG...INFO] number of Q = {i_q_matrix.shape[0]}, P({lambda_i} = {p_vec[i_wl]}, '
              f'S({lambda_i}) = {s_vec[i_wl]}')
        print(f'[Q]: {i_q_matrix[:, 0]}')

        # normalize intensity by K^{lambda} * I^{lambda}(q)
        corrected_intensity_vec = i_q_matrix[:, 2] * k_vec[i_wl]

        # normalize error
        corrected_error_vec = np.zeros_like(corrected_intensity_vec)
        for i_q in range(i_q_matrix.shape[0]):
            #
            print(f'delta I({i_q_matrix[i_q, 0]}, {lambda_i}), I = {i_q_matrix[i_q, 2]}')

            # directly pass the nan values
            if np.isnan(i_q_matrix[i_q, 2]):
                corrected_error_vec[i_q] = np.nan
                print(f'\tNaN')
                continue

            # Using original equations to match with Changwoo's result

            # On reference wavelength (inside qmin and qmax) and thus no correction
            if i_wl == 0:
                # FIXME - need an elegant method to determine minimum wavelength (index)
                # replace qmin_index <= i_q <= qmax_index and lambda_i == ref_wl_vec[i_q, 1]:
                corrected_error_vec[i_q] = i_q_matrix[i_q, 3]**2
                print(f'\tNo correction')
                continue

            # Term 1
            if i_q < qmin_index or i_q > qmax_index:
                # t1 = [delta I(q, wl)]**2 * [P(wl) / S(wl)]**2
                t1_sum = (i_q_matrix[i_q, 3] * p_vec[i_wl] / s_vec[i_wl])**2
                inside = False
            else:
                t1_sum = 0.
                inside = True

            t2_sum = t3_sum = 0.
            print(f'Q     Y     T2     T3')
            for j_q in range(qmin_index, qmax_index + 1):

                # calculate Y
                # Y(q, q', wl) = I(q, wl) * I (q', ref_wl) * S(wl) - I(q, wl) * 2 * I(q', wl) * P(wl)
                y_value = i_q_matrix[i_q, 2] * ref_wl_vec[j_q, 2] * s_vec[i_wl] - i_q_matrix[i_q, 2] * 2. * i_q_matrix[j_q, 2] * p_vec[i_wl]

                if inside and j_q == i_q:
                    # inc = delta I_j(wl_i) * (P^2 * S^2 + 2 * P * S * Y_{q, q}) / S4
                    t1_inc = i_q_matrix[j_q, 3]**2 * (p_vec[i_wl]**2 * s_vec[i_wl]**2 + 2 * p_vec[i_wl] * s_vec[i_wl] * y_value) / s_vec[i_wl]**4
                    t1_sum += t1_inc
                    print(f'[{j_q}]  t1 inc = {t1_inc}')

                # calculate t2_i
                # t2 += [delta I(q', wl)]**2 * Y(q, q'', wl)**2 / S(lw)**4
                t2_inc = i_q_matrix[j_q, 3]**2 * y_value**2 / s_vec[i_wl]**4
                t2_sum += t2_inc

                # calculate t3_i
                # t3: increment = [delta I(q_j, ref_wl[q_j]]^2 * [I(q_j, wl) * I(q, wl)]^2 / S(wl)^2
                # reference: i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 2
                t3_inc = ref_wl_vec[j_q, 3]**2 * i_q_matrix[j_q, 2]**2 * i_q_matrix[i_q, 2]**2 / s_vec[i_wl]**2
                t3_sum += t3_inc

                # DEBUG OUTPUT 1: print(f'{i_q_matrix[i_q, 0]}    {y_value}    {t2_inc}    {t3_inc}')

            # t1 = i_q_matrix[i_q, 3] * p_vec[i_wl] / s_vec[i_wl]
            # print(f'\tterm1 = {t1}')
            #
            # # t2 = sum()
            # t2 = 0.
            # # t3 = sum()
            # t3 = 0.
            # for j_q in range(qmin_index, qmax_index + 1):
            #     # TODO - it is better to restrict the range of j_q to q_min and q_max ...
            #
            #     # t2: increment = delta I(q_j, wl)^2 * [I(q_j, ref_wl(q_j) * S(wl) - 2 * I(q_j, wl) * P(wl)]^2
            #     t2_inc = i_q_matrix[j_q, 3]**2 * (ref_wl_vec[j_q][2] * s_vec[i_wl]
            #                                       - 2 * i_q_matrix[j_q, 2] * p_vec[i_wl])**2
            #     t2 += t2_inc
            #

            #
            #     print(f'\t{j_q}: t2 += {t2_inc}, t3 += {t3_inc}')
            #
            # # error = t1^2 + I(q, wl)^2/S(wl)^4*t2 + I(q, wl)^2/S(wl)^4*t3
            # corrected_error_vec[i_q] = \
            #     t1**2 + i_q_matrix[i_q, 2]**2 / s_vec[i_wl]**4 * t2 + i_q_matrix[i_q, 2]**2 / s_vec[i_wl]**2 * t3

            print(f'WL = {lambda_i}  Q = {i_q_matrix[i_q, 0]}. e^2 = {t1_sum + t2_sum + t3_sum}: '
                  f't1 = {t1_sum}, t2 = {t2_sum}, t3 = {t3_sum}')
            corrected_error_vec[i_q] = t1_sum + t2_sum + t3_sum
            # print(f'sum = {t1_sum + t2_sum + t3_sum}')

            # print(f'error = {corrected_error_vec[i_q]}')
            # raise RuntimeError('DEBUG STOP')

        # update
        new_mod_q = np.concatenate((new_mod_q, i_q_matrix[:, 0]))
        new_q_error = np.concatenate((new_q_error, i_q_matrix[:, 4]))
        new_wavelength = np.concatenate((new_wavelength, i_q_matrix[:, 1]))
        new_intensity = np.concatenate((new_intensity, corrected_intensity_vec))
        new_error_sq = np.concatenate((new_error_sq, corrected_error_vec))

    # construct a new instance of I(Q)
    if delta_q_none:
        new_q_error = None
    corrected_i_of_q = IQmod(new_intensity, np.sqrt(new_error_sq), new_mod_q, new_q_error, new_wavelength)

    return corrected_i_of_q


def correct_intensity_error():
    # correct intensity error: i.e., calcuting new error
    # for a single I(Q, wl)

    # This is based on cleaned up equations optimized for computation efficiency

    # t1 = delta I(q, lambda) * P(lambda) / S(lambda)
    t1 = i_q_matrix[i_q, 3] * p_vec[i_wl] / s_vec[i_wl]
    print(f'\tterm1 = {t1}')

    # t2 = sum()
    t2 = 0.
    # t3 = sum()
    t3 = 0.
    for j_q in range(qmin_index, qmax_index + 1):
        # TODO - it is better to restrict the range of j_q to q_min and q_max ...

        # t2: increment = delta I(q_j, wl)^2 * [I(q_j, ref_wl(q_j) * S(wl) - 2 * I(q_j, wl) * P(wl)]^2
        t2_inc = i_q_matrix[j_q, 3] ** 2 * (ref_wl_vec[j_q][2] * s_vec[i_wl]
                                            - 2 * i_q_matrix[j_q, 2] * p_vec[i_wl]) ** 2
        t2 += t2_inc

        # t3: increment = delta I(q_j, ref_wl[q_j])^2 * I(q_j, wl)^2
        t3_inc = ref_wl_vec[j_q, 3] ** 2 * i_q_matrix[j_q, 2] ** 2
        t3 += t3_inc

        print(f'\t{j_q}: t2 += {t2_inc}, t3 += {t3_inc}')

    # error = t1^2 + I(q, wl)^2/S(wl)^4*t2 + I(q, wl)^2/S(wl)^4*t3
    corrected_error_vec[i_q] = \
        t1 ** 2 + i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 4 * t2 + i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 2 * t3

    print(f't1 = {t1 ** 2}')
    print(f't2 = {i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 4 * t2}')
    print(f't3 = {i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 2 * t3}')

    t1 = t1 ** 2
    t2 = i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 4 * t2
    t3 = i_q_matrix[i_q, 2] ** 2 / s_vec[i_wl] ** 2 * t3
    summed = t1 + t2 + t3
    print(f'sum = {t1 + t2 + t3}')

    return summed
