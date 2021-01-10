# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
# Step 3: Correct wave-dependent incoherence intensity for I(Q, wavelength)
from drtsans.tof.eqsans.elastic_reference_normalization import determine_reference_wavelength_q1d_mesh
import numpy as np


def correct_incoherence_inelastic_1d(i_of_q):
    # This is the envelop method for the complete workflow to correct I(Q1D) counting
    # wavelength-dependent incoherent inelastic

    # Verify: q-bins are same for all wavelength bins: meshgrid between Q and wavelength
    # TODO: based on the requirement as q-bins must be same for all wavelength bins,
    #       it can be pushed backward to the start of the workflow that Q and wavelength are in mesh grid

    # Determine q_min and q_max  that exist in all wavelength for the fitting (minimization) process.
    min_q, max_q = determine_common_mod_q_range(i_of_q)

    # Calculate reference wavelength for each Q: shortest wavelength bin
    ref_wl_matrix = determine_reference_wavelength_q1d(iqmod)

    # Calculate inelastic incoherent factor: b(lambda)
    b_factor_vec = calculate_inelastic_incoherent_factor(i_of_q, min_q, max_q)

    # Optionally select minimum incoherence ...
    if select_min_incoh:
        # re-select minimum incoherence wavelength
        ref_wl_matrix2 = select_min_incoherence(b_factor_vec, unique_wl)
        # re-alculate inelastic incoherent factor: b(lambda)
        b_factor_vec = calculate_inelastic_incoherent_factor(i_of_q, min_q, max_q, ref_wl_matrix2)
        # check: all b_factor_vec[i_wl] shall be greater than ZERO
        assert b_factor_vec > 0

    # Update data for I(Q, wavelength) and delta I(Q, wavelength)


def calculate_b_factors(wl_vec, q_vec, intensity_array, error_array,
                        select_min_incoherence, qmin_index, qmax_index):
    """Determine reference wavelength and then calculate B factor, B error factor.
    
    With option select_min_incoherence, reference wavelength will be reselected according to first
    round of calculation of B
    
    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        wavelength vector
    q_vec: ~numpy.ndarray
        Q vector
    intensity_array: ~numpy.ndarray
        intenisty 2D array
    error_array: ~numpy.ndarray
        intensity error 2D array
    select_min_incoherence: bool
        flag to apply select minimum incoherence algorithm
    qmin_index: int
        index of minimum common q in q vector
    qmax_index: int
        index of maximum common q in q vector (included)

    Returns
    -------

    """
    # Sanity check
    assert intensity_array.shape == error_array.shape
    assert wl_vec.shape[0] == intensity_array.shape[1]
    assert q_vec.shape[0] == error_array.shape[0]

    # determine the reference wavelength: minimum wavelength bin
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, intensity_array, error_array,
                                                        qmin_index, qmax_index, 0)

    # calculate b(lambda_i) and delta b(lambda_i) if it is final
    b_array = calculate_b_error_b(wl_vec, intensity_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                  calculate_b_error=not select_min_incoherence)
    print(f'[DEVELOP OUTPUT] First round B vector: {b_array[0]}')

    # If JSON parameter “selectMinIncoh” is true
    if select_min_incoherence and np.argmin(b_array[0]) > 0:
        # reselect reference wavelength to the wavelength bin with minimum b
        ref_wl_index = np.argmin(b_array[0])
        print(f'[DEBUG...] re-selected reference wavelength index = {ref_wl_index}')
        # (re)determine the reference wavelengths' intensities and errors
        ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, intensity_array, error_array,
                                                            qmin_index, qmax_index, ref_wl_index)
        # (re)calculate b array
        b_array = calculate_b_error_b(wl_vec, intensity_array, error_array, qmin_index, qmax_index, ref_wl_ie,
                                      calculate_b_error=True)
        # verify
        assert b_array[np.isfinite(b_array)].min() >= -1E-20, f'B array has negative values: {b_array}'

    return b_array, ref_wl_ie


def calculate_b_error_b(wl_vec, intensity_array, error_array, qmin_index, qmax_index,
                        ref_wavelengths, calculate_b_error):
    """

    Parameters
    ----------
    wl_vec
    q_vec
    intensity_array
    error_array
    qmin_index
    qmax_index: int
        index of maximum common Q (included)
    ref_wavelengths
    calculate_b_error

    Returns
    -------
    ~numpy.ndarray
        row 0: B factor, row 1: delta B

    """
    # Declare B factor array
    b_factor_array = np.zeros(shape=(2, len(wl_vec)), dtype='float')

    # Calculate B factors
    # b[wl] = - 1/N sum_{q_k=q_min}^{q_max} [RefI(q_k) - I(q_k, wl)]
    num_q = qmax_index + 1 - qmin_index
    # operation into a (num_q, num_wl) 2D array
    b_vec = ref_wavelengths.intensity_vec[qmin_index:qmax_index+1].reshape((num_q, 1))\
            - intensity_array[qmin_index:qmax_index+1, :]
    b_factor_array[0] = -1. / num_q * np.sum(b_vec, axis=0)

    # Calculate B error (delta B) as an option
    if calculate_b_error:
        # delta b(wl)^2 = 1/N^2 sum_{q_k=qmin}^{qmax} [(delta I(q_k, ref_wl))^2 + (delta I(q_k, wl))^2]
        # operation into a (num_q, num_wl) 2D array
        b2_vec = (ref_wavelengths.error_vec[qmin_index:qmax_index+1].reshape((num_q, 1)))**2\
                 + (error_array[qmin_index:qmax_index+1, :])**2
        b_factor_array[1] = 1. / num_q * np.sqrt(b2_vec.sum(axis=0))

    return b_factor_array


def correct_intensity_error(wavelength_vec, q_vec, intensity_array, error_array, b_vector,
                            qmin_index, qmax_index, ref_wl_ie):

    # Sanity checks
    assert intensity_array.shape == error_array.shape
    assert wavelength_vec.shape[0] == intensity_array.shape[1]
    assert q_vec.shape[0] == error_array.shape[0]
    assert b_vector.shape[1] == wavelength_vec.shape[0]
    num_q = q_vec.shape[0]

    # Correct intensity: I'(q, wl_i) = I(q, wl_i) - b(wl_i)
    corrected_intensities = intensity_array[:, :] - b_vector[0].reshape([num_q, 1])

    # Correct error:
    num_common_q = qmax_index - qmin_index + 1
    corrected_errors = np.zeros_like(error_array)

    # Loop over wavelength
    for i_wl in range(wavelength_vec.shape[0]):

        # outside q_min and q_max
        # term1[q_j] = (error[q_j, wl]^2
        term1 = error_array[:, i_wl]**2

        # term2 = 1/N^2 sum_{q_k=q_min:qmax}[RefError(q_k)^2 + error(q_k, wl)^2]
        # term2 is a single value
        term2 = 1. / num_common_q**2 * np.sum(ref_wl_ie.error_vec[qmin_index:qmax_index+1]**2 +
                                              error_array[qmin_index:qmax_index+1]**2)
        print(f'Term2 = {term2}')

        # make correct term1 for those inside qmin and qmax
        # term1[q_j] = (error[q_j, wl]^2 * (1 - 2/N)
        term1[qmin_index:qmax_index+1] *= (1 - 2/num_common_q)

        # sum
        term1 += term2
        print(f'Term1 shape = {term1.shape}')

        # assign to output
        corrected_errors[:, i_wl] = term1.reshape((num_q, 1))

    return corrected_intensities, np.sqrt(corrected_errors)
