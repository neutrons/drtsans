# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from drtsans.dataobjects import IQazimuthal, IQmod, save_iqmod, verify_same_q_bins

__all__ = [
    "normalize_by_elastic_reference_all",
    "normalize_by_elastic_reference_1d",
    "normalize_by_elastic_reference_2d",
    "determine_reference_wavelength_q1d_mesh",
    "reshape_q_wavelength_matrix",
    "build_i_of_q1d",
    "determine_common_mod_q_range_mesh",
]


@dataclass
class ReferenceWavelengths:
    """
    Class for keeping track of reference wavelength for each momentum transfer Q (1D)

    Parameters
    ----------
    q_values: ~numpy.ndarray
        vector for Q
    ref_wavelengths: ~numpy.ndarray
        vector for reference wavelength vector for each Q
    intensities: ~numpy.ndarray
        vector for intensities of (Q, reference wavelength)
    errors: ~numpy.ndarray
        vector for errors of (Q, reference wavelength)
    """

    q_vec: np.ndarray
    ref_wl_vec: np.ndarray
    intensity_vec: np.ndarray
    error_vec: np.ndarray


def reshape_q_wavelength_matrix(i_of_q: IQmod) -> tuple:
    """Reshape I(Q) into a mesh grid of (Q, wavelength)

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from

    Returns
    -------
    tuple
        wavelength vector, q vector,  intensity (2D), error (2D), dq array (2D) or None

    """
    # Create a matrix for q, wavelength, intensity and error
    if i_of_q.delta_mod_q is None:
        i_q_wl_matrix = np.array([i_of_q.mod_q, i_of_q.wavelength, i_of_q.intensity, i_of_q.error])
    else:
        i_q_wl_matrix = np.array(
            [
                i_of_q.mod_q,
                i_of_q.wavelength,
                i_of_q.intensity,
                i_of_q.error,
                i_of_q.delta_mod_q,
            ]
        )
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
    if i_of_q.delta_mod_q is not None:
        dq_array = i_q_wl_matrix[:, 4].reshape((q_vector.shape[0], wl_vector.shape[0]))
    else:
        dq_array = None

    return wl_vector, q_vector, intensity_array, error_array, dq_array


def normalize_by_elastic_reference_all(
    i_of_q_2d, i_of_q_1d, ref_i_of_q_1d, output_wavelength_dependent_profile=False, output_dir=None
):
    """Normalize I(Q2D) and I(Q1D) by elastic reference run

    Parameters
    ----------
    i_of_q_2d: ~drtsans.dataobjects.IQazimuthal
        Input I(Q2D, wavelength) to normalize
    i_of_q_1d: ~drtsans.dataobjects.IQmod
        Input I(Q1D, wavelength) to normalize
    ref_i_of_q_1d: ~drtsans.dataobjects.IQmod
        Input I(Q1D, wavelength) as elastic reference run
    output_wavelength_dependent_profile: bool
        If True then output Iq for each wavelength before and after k correction
    output_dir: str
        output directory for Iq profiles

    Returns
    -------
    tuple
        normalized I(Q2D), normalized I(Q1D), K vector and delta K vector
    """
    # check i_of_q and ref_i_of_q shall have same binning
    if not verify_same_q_bins(i_of_q_1d, ref_i_of_q_1d, False, tolerance=1e-3):
        raise RuntimeError("Input I(Q) and elastic reference I(Q) have different Q and wavelength binning")

    wl_vec = np.unique(i_of_q_1d.wavelength)
    q_vec = np.unique(i_of_q_1d.mod_q)

    # Calculate the normalization factor for each wavelength
    k_vec, k_error_vec = calculate_elastic_reference_normalization(wl_vec, q_vec, ref_i_of_q_1d)

    # 1D normalization
    iq1d_wl = normalize_by_elastic_reference_1d(
        i_of_q_1d,
        k_vec,
        k_error_vec,
        output_wavelength_dependent_profile,
        output_dir,
    )

    # 2D normalization
    iq2d_wl = normalize_by_elastic_reference_2d(i_of_q_2d, k_vec, k_error_vec)

    return iq2d_wl, iq1d_wl, k_vec, k_error_vec


def calculate_elastic_reference_normalization(wl_vec, q_vec, ref_i_of_q):
    """Calculate the elastic reference normalization factor (K) for each wavelength

    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        Vector of wavelengths in I(Q) to normalize
    q_vec: ~numpy.ndarray
        Vector of Q:s in I(Q) to normalize
    ref_i_of_q: ~drtsans.dataobjects.IQmod
        Elastic reference run I(Q, wavelength)

    Returns
    -------
    tuple
        K vector and delta K vector

    """
    # Reshape Q, wavelength, intensities and errors to unique 1D array or 2D array
    _, _, ref_i_array, ref_error_array, _ = reshape_q_wavelength_matrix(ref_i_of_q)

    # Calculate Qmin and Qmax
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, ref_i_array)

    # Calculate reference
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(
        wl_vec, q_vec, ref_i_array, ref_error_array, qmin_index, qmax_index
    )

    # Calculate scale factor
    k_vec, k_error_vec = calculate_scale_factor_mesh_grid(
        wl_vec, ref_i_array, ref_error_array, ref_wl_ie, qmin_index, qmax_index
    )

    return k_vec, k_error_vec


def normalize_by_elastic_reference_1d(
    i_of_q: IQmod,
    k_vec: np.ndarray,
    k_error_vec: np.ndarray,
    output_wavelength_dependent_profile: bool = False,
    output_dir: Optional[str] = None,
) -> IQmod:
    """Normalize I(Q1D) by wavelength-dependent elastic reference normalization factor

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to normalize
    k_vec: ~numpy.ndarray
        Elastic reference normalization factors (one for each wavelength)
    k_error_vec: ~numpy.ndarray
        Elastic reference normalization factor errors (one for each wavelength)
    output_wavelength_dependent_profile: bool
        If True then output Iq for each wavelength before and after k correction
    output_dir: str
        output directory for Iq profiles

    Returns
    -------
    tuple
        normalized Q(1D), K vector and delta K vector

    """
    # Reshape Q, wavelength, intensities and errors to unique 1D array or 2D array
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_of_q)

    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_before_k_correction.dat")
            save_iqmod(
                IQmod(
                    intensity=i_array[:, tmpwlii],
                    error=error_array[:, tmpwlii],
                    mod_q=q_vec,
                    delta_mod_q=dq_array[:, tmpwlii],
                ),
                tmpfn,
            )

    # Normalize
    normalized = normalize_intensity_q1d(
        wl_vec,
        q_vec,
        i_array,
        error_array,
        k_vec,
        k_error_vec,
    )

    # Convert normalized intensities and errors to IModQ
    normalized_i_of_q = build_i_of_q1d(wl_vec, q_vec, normalized[0], normalized[1], dq_array)

    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_after_k_correction.dat")
            save_iqmod(
                IQmod(
                    intensity=normalized[0][:, tmpwlii],
                    error=normalized[1][:, tmpwlii],
                    mod_q=q_vec,
                    delta_mod_q=dq_array[:, tmpwlii],
                ),
                tmpfn,
            )

    return normalized_i_of_q


def normalize_by_elastic_reference_2d(i_of_q, k_vec, k_error_vec):
    """Normalize I(Q2D) by wavelength-dependent elastic reference normalization factor

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        Input I(Q2D, wavelength) to normalize
    k_vec: ~numpy.ndarray
        Elastic reference normalization factors (one for each wavelength)
    k_error_vec: ~numpy.ndarray
        Elastic reference normalization factor errors (one for each wavelength)

    Returns
    -------
    ~drtsans.dataobjects.IQazimuthal
        normalized I(Q2D)
    """
    intensity_array = i_of_q.intensity
    error_array = i_of_q.error

    # Reshape vectors to be easily indexed by wavelength
    num_wl = len(np.unique(i_of_q.wavelength))
    sizeX = i_of_q.qx.shape[0]
    sizeY = i_of_q.qy.shape[0]
    intensity_3d = intensity_array.transpose().reshape((num_wl, sizeX, sizeY))
    error_3d = error_array.transpose().reshape((num_wl, sizeX, sizeY))

    # Scale each wavelength by the corresponding normalization factor, K
    normalized_intensity_array = np.zeros_like(intensity_3d)
    for i_wl in range(num_wl):
        normalized_intensity_array[i_wl] = intensity_3d[i_wl] * k_vec[i_wl]

    # Lowest wavelength bin does not require normalization as K = 1, i_wl = 0
    normalized_error2_array = np.zeros_like(error_3d)
    normalized_error2_array[0, :, :] = error_3d[0, :, :] ** 2

    for i_wl in range(1, num_wl):
        # Calculate error as dI^2 = dK^2 * I^2 + K^2 * dI^2
        normalized_error2_array[i_wl] = (
            k_error_vec[i_wl] ** 2 * intensity_3d[i_wl] ** 2 + k_vec[i_wl] ** 2 * error_3d[i_wl] ** 2
        )
    normalized_error_array = np.sqrt(normalized_error2_array)

    # Transform vectors back to their original shape
    normalized_intensity_array = normalized_intensity_array.reshape((sizeX * num_wl, sizeY)).transpose()
    normalized_error_array = normalized_error_array.reshape((sizeX * num_wl, sizeY)).transpose()

    normalized_i_of_q = IQazimuthal(
        intensity=normalized_intensity_array,
        error=normalized_error_array,
        qx=i_of_q.qx,
        qy=i_of_q.qy,
        wavelength=i_of_q.wavelength,
        delta_qx=i_of_q.delta_qx,
        delta_qy=i_of_q.delta_qy,
    )

    return normalized_i_of_q


def build_i_of_q1d(wl_vector, q_vector, intensity_array, error_array, delta_q_array) -> IQmod:
    """From wavelength, Q, intensity, error and delta Q to build I(Q1D)

    This is the reversed operation to method reshape_q_wavelength_matrix

    Parameters
    ----------
    wl_vector: ~numpy.ndarray
        wavelength (1D)
    q_vector: ~numpy.ndarray
        Q (1D)
    intensity_array: ~numpy.ndarray
        intensities (2D)
    error_array: ~numpy.ndarray
        intensity errors (2D)
    delta_q_array: ~numpy.ndarray
        delta Q (1D) size = number wavelength * number Q

    Returns
    -------
    ~drtsans.dataobjects.IQmod
        constructed I(Q, wavelength)

    """
    # assume that intensity, error and delta q have the same as (num_q, num_wl)
    assert intensity_array.shape[0] == q_vector.shape[0] and intensity_array.shape[1] == wl_vector.shape[0]

    # tile wave length
    wl_array_1d = np.tile(wl_vector, q_vector.shape[0])
    q_array_1d = np.repeat(q_vector, wl_vector.shape[0])

    # flatten intensity, error and optionally delta q
    intensity_array = intensity_array.flatten()
    error_array = error_array.flatten()
    if delta_q_array is not None:
        delta_q_array = delta_q_array.flatten()

    return IQmod(
        intensity=intensity_array,
        error=error_array,
        mod_q=q_array_1d,
        wavelength=wl_array_1d,
        delta_mod_q=delta_q_array,
    )


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
    assert q_vec.shape[0] == intensity_array.shape[0], "Shape mismatch"

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
        raise RuntimeError("Unable to find common q range")

    return qmin_index, qmax_index


def calculate_scale_factor_mesh_grid(wl_vec, intensity_array, error_array, ref_wl_intensities, qmin_index, qmax_index):
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
        K vector, K error vector
    """
    # Check input
    assert wl_vec.shape[0] == intensity_array.shape[1]

    k_vec = np.zeros_like(wl_vec)
    k_error2_vec = np.zeros_like(wl_vec)

    for i_wl, lambda_i in enumerate(wl_vec):
        # P(wl) = sum_q I(q, ref_wl) * I(q, wl)
        p_value = np.sum(
            ref_wl_intensities.intensity_vec[qmin_index : qmax_index + 1]
            * intensity_array[:, i_wl][qmin_index : qmax_index + 1]
        )
        # S(wl) = sum_q I(q, wl)**2
        s_value = np.sum(intensity_array[:, i_wl][qmin_index : qmax_index + 1] ** 2)

        # Calculate K(wl)
        k_vec[i_wl] = p_value / s_value

        # Calculate K(wl) error
        term0 = error_array[:, i_wl][qmin_index : qmax_index + 1]
        term1 = (
            ref_wl_intensities.intensity_vec[qmin_index : qmax_index + 1] * s_value
            - 2.0 * intensity_array[:, i_wl][qmin_index : qmax_index + 1] * p_value
        ) / s_value**2
        term2 = ref_wl_intensities.error_vec[qmin_index : qmax_index + 1]
        term3 = intensity_array[:, i_wl][qmin_index : qmax_index + 1] / s_value

        k_error2_vec[i_wl] = np.sum((term0 * term1) ** 2 + (term2 * term3) ** 2)

    return k_vec, np.sqrt(k_error2_vec)


def determine_reference_wavelength_q1d_mesh(
    wavelength_vec,
    q_vec,
    intensity_array,
    error_array,
    qmin_index,
    qmax_index,
    min_wl_index=0,
):
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
    assert wavelength_vec.shape[0] == intensity_array.shape[1], (
        f"Wavelength dimension = {wavelength_vec.shape} ," f"Intensity  dimension = {intensity_array.shape}"
    )

    # Minimum wavelength bin is the reference wavelength
    min_wl_vec = np.zeros_like(q_vec) + wavelength_vec[min_wl_index]

    # Minimum intensity and error
    min_intensity_vec = np.copy(intensity_array[:, min_wl_index])
    min_error_vec = np.copy(error_array[:, min_wl_index])

    # Set the unused defined reference wavelength (outside of qmin and qmax)'s
    # intensity and error to nan
    min_intensity_vec[0:qmin_index] = np.nan
    min_intensity_vec[qmax_index + 1 :] = np.nan
    min_error_vec[0:qmin_index] = np.nan
    min_error_vec[qmax_index + 1 :] = np.nan

    return ReferenceWavelengths(q_vec, min_wl_vec, min_intensity_vec, min_error_vec)


def normalize_intensity_q1d(
    wl_vec,
    q_vec,
    intensity_array,
    error_array,
    k_vec,
    k_error_vec,
):
    """Normalize Q1D intensities and errors

    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        1D vector of wavelength (in ascending order)
    q_vec: ~numpy.ndarray
        1D vector of Q (in ascending order)
    intensity_array: ~numpy.ndarray
        2D array of intensities, shape[0] = number of Q, shape[1] = number of wavelength
    error_array: ~numpy.ndarray
        2D array of errors, shape[0] = number of Q, shape[1] = number of wavelength
    k_vec: ~numpy.ndarray
        calculated K vector
    k_error_vec: ~numpy.ndarray
        calculated K error vector

    Returns
    -------
    tuple
        normalized I(Q1D), normalized error(Q1D)

    """

    # Sanity check
    assert wl_vec.shape[0] == intensity_array.shape[1]  # wavelength as lambda
    assert q_vec.shape[0] == error_array.shape[0]  # points as q
    assert intensity_array.shape == error_array.shape

    # Normalized intensities
    normalized_intensity_array = intensity_array * k_vec
    normalized_error2_array = np.zeros_like(error_array)

    # Lowest wavelength bin does not require normalization as K = 1, i_wl = 0
    normalized_error2_array[:, 0] = error_array[:, 0] ** 2

    # Loop over wavelength
    num_wl = wl_vec.shape[0]
    for i_wl in range(1, num_wl):
        # Calculate error as dI^2 = dK^2 * I^2 + K^2 * dI^2
        normalized_error2_array[:, i_wl] = (
            k_error_vec[i_wl] ** 2 * intensity_array[:, i_wl] ** 2 + k_vec[i_wl] ** 2 * error_array[:, i_wl] ** 2
        )

    return normalized_intensity_array, np.sqrt(normalized_error2_array)
