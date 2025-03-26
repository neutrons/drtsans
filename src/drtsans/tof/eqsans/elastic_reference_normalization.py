# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from drtsans.dataobjects import IQazimuthal, IQmod, verify_same_q_bins, I1DAnnular, save_i1d, getDataType, DataType

__all__ = [
    "normalize_by_elastic_reference_all",
    "normalize_by_elastic_reference_1d",
    "normalize_by_elastic_reference_2d",
    "determine_reference_wavelength_intensity_mesh",
    "reshape_intensity_domain_meshgrid",
    "build_i1d_from_intensity_domain_meshgrid",
    "build_i1d_one_wl_from_intensity_domain_meshgrid",
    "determine_common_domain_range_mesh",
]


@dataclass
class ReferenceWavelengths:
    """
    Class for keeping track of reference wavelength for each bin (Q or phi)

    Parameters
    ----------
    x_vec: ~numpy.ndarray
        vector for Q or phi
    ref_wl_vec: ~numpy.ndarray
        vector for reference wavelength vector for each Q or phi
    intensity_vec: ~numpy.ndarray
        vector for intensities of (Q, reference wavelength) or (phi, reference wavelength)
    error_vec: ~numpy.ndarray
        vector for errors of (Q, reference wavelength) or (phi, reference wavelength)
    """

    x_vec: np.ndarray
    ref_wl_vec: np.ndarray
    intensity_vec: np.ndarray
    error_vec: np.ndarray


def reshape_intensity_domain_meshgrid(i1d: IQmod | I1DAnnular) -> tuple:
    """Reshape I(Q)/I(phi) into a mesh grid of (Q, wavelength) or (phi, wavelength)

    Parameters
    ----------
    i1d: IQmod | I1DAnnular
        Input I(Q, wavelength) or I(phi, wavelength) to find common range from

    Returns
    -------
    tuple
        wavelength vector, Q/phi vector, intensity (2D), error (2D), dq array (2D) or None

    """
    # Create a matrix for Q/phi, wavelength, intensity and error
    # x is either momentum transfer, Q, or azimuthal angle, phi
    if i1d.delta_x is None:
        i_x_wl_matrix = np.array([i1d.x, i1d.wavelength, i1d.intensity, i1d.error])
    else:
        i_x_wl_matrix = np.array(
            [
                i1d.x,
                i1d.wavelength,
                i1d.intensity,
                i1d.error,
                i1d.delta_x,
            ]
        )
    i_x_wl_matrix = i_x_wl_matrix.transpose()

    # Order by wavelength and x
    i_x_wl_matrix = i_x_wl_matrix[np.lexsort((i_x_wl_matrix[:, 1], i_x_wl_matrix[:, 0]))]

    # Unique wavelength and unique x
    wl_vector = np.unique(i1d.wavelength)
    x_vector = np.unique(i1d.x)
    # verify whether (x, wl) are on mesh grid by checking unique x and wavelength
    assert wl_vector.shape[0] * x_vector.shape[0] == i1d.intensity.shape[0]

    # Reformat
    intensity_array = i_x_wl_matrix[:, 2].reshape((x_vector.shape[0], wl_vector.shape[0]))
    error_array = i_x_wl_matrix[:, 3].reshape((x_vector.shape[0], wl_vector.shape[0]))
    if i1d.delta_x is not None:
        delta_x_array = i_x_wl_matrix[:, 4].reshape((x_vector.shape[0], wl_vector.shape[0]))
    else:
        delta_x_array = None

    return wl_vector, x_vector, intensity_array, error_array, delta_x_array


def normalize_by_elastic_reference_all(
    i_of_q_2d, i_1d, ref_i_1d, output_wavelength_dependent_profile=False, output_dir=None
):
    """Normalize I(Q2D) and I(1D) by elastic reference run

    Parameters
    ----------
    i_of_q_2d: ~drtsans.dataobjects.IQazimuthal
        Input I(Q2D, wavelength) to normalize
    i_1d: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Input I(Q1D, wavelength) or I(phi, wavelength) to normalize
    ref_i_1d: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Input I(Q1D, wavelength) or I(phi, wavelength) as elastic reference run
    output_wavelength_dependent_profile: bool
        If True then output I(1D) for each wavelength before and after k correction
    output_dir: str
        output directory for I(1D) profiles

    Returns
    -------
    tuple
        normalized I(Q2D), normalized I(1D), K vector and delta K vector
    """
    # check i_of_q and ref_i_of_q shall have same binning
    if not verify_same_q_bins(i_1d, ref_i_1d, False, tolerance=1e-3):
        raise RuntimeError("Input I(1D) and elastic reference I(1D) have different binning")

    wl_vec = np.unique(i_1d.wavelength)
    # x is either momentum transfer, Q, or azimuthal angle, phi
    x_vec = np.unique(i_1d.x)

    # Calculate the normalization factor for each wavelength
    k_vec, k_error_vec = calculate_elastic_reference_normalization(wl_vec, x_vec, ref_i_1d)

    # 1D normalization
    i1d_wl = normalize_by_elastic_reference_1d(
        i_1d,
        k_vec,
        k_error_vec,
        output_wavelength_dependent_profile,
        output_dir,
    )

    # 2D normalization
    iq2d_wl = normalize_by_elastic_reference_2d(i_of_q_2d, k_vec, k_error_vec)

    return iq2d_wl, i1d_wl, k_vec, k_error_vec


def calculate_elastic_reference_normalization(wl_vec, x_vec, ref_i_1d):
    """Calculate the elastic reference normalization factor (K) for each wavelength

    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        Vector of wavelengths in I(1D) to normalize
    x_vec: ~numpy.ndarray
        Vector of Q:s, or phi:s, in I(1D) to normalize
    ref_i_1d: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Elastic reference run I(Q, wavelength) or I(phi, wavelength)

    Returns
    -------
    tuple
        K vector and delta K vector

    """
    # Reshape Q/phi, wavelength, intensities and errors to unique 1D array or 2D array
    _, _, ref_i_array, ref_error_array, _ = reshape_intensity_domain_meshgrid(ref_i_1d)

    # Calculate Qmin and Qmax
    qmin_index, qmax_index = determine_common_domain_range_mesh(x_vec, ref_i_array)

    # Calculate reference
    ref_wl_ie = determine_reference_wavelength_intensity_mesh(
        wl_vec, x_vec, ref_i_array, ref_error_array, qmin_index, qmax_index
    )

    # Calculate scale factor
    k_vec, k_error_vec = calculate_scale_factor_mesh_grid(
        wl_vec, ref_i_array, ref_error_array, ref_wl_ie, qmin_index, qmax_index
    )

    return k_vec, k_error_vec


def normalize_by_elastic_reference_1d(
    i1d: IQmod | I1DAnnular,
    k_vec: np.ndarray,
    k_error_vec: np.ndarray,
    output_wavelength_dependent_profile: bool = False,
    output_dir: Optional[str] = None,
) -> IQmod | I1DAnnular:
    """Normalize I(Q1D) or I(phi) by wavelength-dependent elastic reference normalization factor

    Parameters
    ----------
    i1d: IQmod | I1DAnnular
        Input I(Q, wavelength) or I(phi, wavelength) to normalize
    k_vec: ~numpy.ndarray
        Elastic reference normalization factors (one for each wavelength)
    k_error_vec: ~numpy.ndarray
        Elastic reference normalization factor errors (one for each wavelength)
    output_wavelength_dependent_profile: bool
        If True then output I for each wavelength before and after k correction
    output_dir: str
        output directory for intensity profiles

    Returns
    -------
    tuple
        normalized I(Q1D) or I(phi), K vector and delta K vector

    """
    # Reshape Q/phi, wavelength, intensities and errors to unique 1D array or 2D array
    wl_vec, x_vec, i_array, error_array, dq_array = reshape_intensity_domain_meshgrid(i1d)

    i1d_type = getDataType(i1d)
    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_before_k_correction.dat")
            i1d_wl = build_i1d_one_wl_from_intensity_domain_meshgrid(
                x_vec, i_array, error_array, dq_array, tmpwlii, i1d_type
            )
            save_i1d(i1d_wl, tmpfn)

    # Normalize
    normalized = normalize_intensity_1d(
        wl_vec,
        x_vec,
        i_array,
        error_array,
        k_vec,
        k_error_vec,
    )

    # Convert normalized intensities and errors to object of type `i1d_type`
    normalized_i1d = build_i1d_from_intensity_domain_meshgrid(
        wl_vec, x_vec, normalized[0], normalized[1], dq_array, i1d_type
    )

    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_after_k_correction.dat")
            i1d_wl = build_i1d_one_wl_from_intensity_domain_meshgrid(
                x_vec, normalized[0], normalized[1], dq_array, tmpwlii, i1d_type
            )
            save_i1d(i1d_wl, tmpfn)

    return normalized_i1d


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


def build_i1d_from_intensity_domain_meshgrid(
    wl_vector, x_vector, intensity_array, error_array, delta_q_array, i1d_type
) -> IQmod | I1DAnnular:
    """Build I(1D) from a mesh grid representation of intensity, error and delta Q
    as a function of Q (or phi) and wavelength

    This is the reversed operation to method
    :py:meth:`~drtsans.tof.eqsans.elastic_reference_normalization.reshape_intensity_domain_meshgrid`

    Parameters
    ----------
    wl_vector: ~numpy.ndarray
        wavelength (1D)
    x_vector: ~numpy.ndarray
        Q or phi (1D)
    intensity_array: ~numpy.ndarray
        intensities (2D)
    error_array: ~numpy.ndarray
        intensity errors (2D)
    delta_q_array: ~numpy.ndarray
        delta Q (1D) size = number wavelength * number Q, not used for i1d_type = DataType.I_ANNULAR
    i1d_type: DataType
        output data type for the intensity profile

    Returns
    -------
    ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        constructed I(Q, wavelength) or I(phi, wavelength)

    """
    # assume that intensity, error and delta q have the same as (num_q, num_wl)
    assert intensity_array.shape[0] == x_vector.shape[0] and intensity_array.shape[1] == wl_vector.shape[0]

    # tile wave length
    wl_array_1d = np.tile(wl_vector, x_vector.shape[0])
    q_array_1d = np.repeat(x_vector, wl_vector.shape[0])

    # flatten intensity, error and optionally delta q
    intensity_array = intensity_array.flatten()
    error_array = error_array.flatten()
    if delta_q_array is not None:
        delta_q_array = delta_q_array.flatten()

    if i1d_type == DataType.IQ_MOD:
        i1d = IQmod(
            intensity=intensity_array,
            error=error_array,
            mod_q=q_array_1d,
            wavelength=wl_array_1d,
            delta_mod_q=delta_q_array,
        )
    elif i1d_type == DataType.I_ANNULAR:
        i1d = I1DAnnular(
            intensity=intensity_array,
            error=error_array,
            phi=q_array_1d,
            wavelength=wl_array_1d,
        )
    else:
        raise TypeError(f"{i1d_type} is not supported")

    return i1d


def build_i1d_one_wl_from_intensity_domain_meshgrid(
    x_vector, intensity_array, error_array, delta_q_array, wl_index, i1d_type
) -> IQmod | I1DAnnular:
    """Build I(1D) for one wavelength from a mesh grid representation of intensity,
    error and delta Q as a function of Q (or phi) and wavelength

    Parameters
    ----------
    x_vector: ~numpy.ndarray
        Q or phi (1D)
    intensity_array: ~numpy.ndarray
        intensities (2D)
    error_array: ~numpy.ndarray
        intensity errors (2D)
    delta_q_array: ~numpy.ndarray
        delta Q (1D) size = number wavelength * number Q, not used for i1d_type = DataType.I_ANNULAR
    wl_index: int
        Wavelength axis index to extract
    i1d_type: DataType
        output data type for the intensity profile

    Returns
    -------
    ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        constructed I(Q, wavelength) or I(phi, wavelength)

    """
    if i1d_type == DataType.IQ_MOD:
        return IQmod(
            intensity=intensity_array[:, wl_index],
            error=error_array[:, wl_index],
            mod_q=x_vector,
            delta_mod_q=delta_q_array[:, wl_index],
        )
    elif i1d_type == DataType.I_ANNULAR:
        return I1DAnnular(
            intensity=intensity_array[:, wl_index],
            error=error_array[:, wl_index],
            phi=x_vector,
        )
    else:
        raise TypeError(f"{i1d_type} is not supported")


def determine_common_domain_range_mesh(x_vec, intensity_array):
    """Determine the common Q or phi range among all the wavelengths such that I(Q/phi, lambda) does exist.

    This method assumes that I(Q/phi, wavelength) are on mesh grid of Q/phi and wavelength

    Detailed requirement:
        Determine q_min and q_max that exist in all I(q, lambda) for the fitting (minimization) process

    Parameters
    ----------
    x_vec: numpy.ndarray
        vector of sorted unique Q or phi, depending on whether the 1D binning is scalar or annular
    intensity_array: numpy.ndarray
        2D array of intensity.  Each row is of same wavelength

    Returns
    -------
    tuple
        index of x_min and x_max

    """
    xmin_index = None
    xmax_index = None

    # Sanity check
    assert x_vec.shape[0] == intensity_array.shape[0], "Shape mismatch"

    num_x = x_vec.shape[0]
    for x_index in range(num_x):
        if len(np.where(np.isnan(intensity_array[x_index]))[0]) == 0:
            xmin_index = x_index
            break
    for x_index in range(num_x - 1, -1, -1):
        if len(np.where(np.isnan(intensity_array[x_index]))[0]) == 0:
            xmax_index = x_index
            break

    if xmin_index is None:
        raise RuntimeError("Unable to find common range")

    return xmin_index, xmax_index


def calculate_scale_factor_mesh_grid(wl_vec, intensity_array, error_array, ref_wl_intensities, xmin_index, xmax_index):
    """Same functionality as calculate_scale_factor but the algorithm is improved
    as I(Q, wavelength) or I(phi, wavelength) are in meshgrid

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
    xmin_index: int
        index of min x in x vector (actually Q or phi)
    xmax_index: int
        index of max x in x vector (actually Q or phi)

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
            ref_wl_intensities.intensity_vec[xmin_index : xmax_index + 1]
            * intensity_array[:, i_wl][xmin_index : xmax_index + 1]
        )
        # S(wl) = sum_q I(q, wl)**2
        s_value = np.sum(intensity_array[:, i_wl][xmin_index : xmax_index + 1] ** 2)

        # Calculate K(wl)
        k_vec[i_wl] = p_value / s_value

        # Calculate K(wl) error
        term0 = error_array[:, i_wl][xmin_index : xmax_index + 1]
        term1 = (
            ref_wl_intensities.intensity_vec[xmin_index : xmax_index + 1] * s_value
            - 2.0 * intensity_array[:, i_wl][xmin_index : xmax_index + 1] * p_value
        ) / s_value**2
        term2 = ref_wl_intensities.error_vec[xmin_index : xmax_index + 1]
        term3 = intensity_array[:, i_wl][xmin_index : xmax_index + 1] / s_value

        k_error2_vec[i_wl] = np.sum((term0 * term1) ** 2 + (term2 * term3) ** 2)

    return k_vec, np.sqrt(k_error2_vec)


def determine_reference_wavelength_intensity_mesh(
    wavelength_vec,
    x_vec,
    intensity_array,
    error_array,
    xmin_index,
    xmax_index,
    min_wl_index=0,
):
    """Determine the reference wavelength for each x (Q or phi).

    The reference wavelength of a specific Q or (qx, qy) or phi
    is defined as the shortest wavelength for all the finite I(Q, wavelength) or
    I(qx, qy, wavelength) or I(phi, wavelength)

    Parameters
    ----------
    wavelength_vec: numpy.ndarray
        Wavelength vector
    x_vec: numpy.ndarray
        Vector of Q or phi
    intensity_array: numpy.ndarray
        Intensity array (2D)
    error_array: numpy.ndarray
        Error array (2D)
    xmin_index: int
        index of xmin in x-vector
    xmax_index: int
        index of xmax in x-vector

    Returns
    -------
    ReferenceWavelengths
        Reference wavelengths for each momentum transfer Q, or azimuthal angle phi, and the
        corresponding intensity and error

    """
    # Sanity check
    assert wavelength_vec.shape[0] == intensity_array.shape[1], (
        f"Wavelength dimension = {wavelength_vec.shape} ,Intensity  dimension = {intensity_array.shape}"
    )

    # Minimum wavelength bin is the reference wavelength
    min_wl_vec = np.zeros_like(x_vec) + wavelength_vec[min_wl_index]

    # Minimum intensity and error
    min_intensity_vec = np.copy(intensity_array[:, min_wl_index])
    min_error_vec = np.copy(error_array[:, min_wl_index])

    # Set the unused defined reference wavelength (outside of xmin and xmax)'s
    # intensity and error to nan
    min_intensity_vec[0:xmin_index] = np.nan
    min_intensity_vec[xmax_index + 1 :] = np.nan
    min_error_vec[0:xmin_index] = np.nan
    min_error_vec[xmax_index + 1 :] = np.nan

    return ReferenceWavelengths(x_vec, min_wl_vec, min_intensity_vec, min_error_vec)


def normalize_intensity_1d(
    wl_vec,
    x_vec,
    intensity_array,
    error_array,
    k_vec,
    k_error_vec,
):
    """Normalize 1D intensities and errors

    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        1D vector of wavelength (in ascending order)
    x_vec: ~numpy.ndarray
        1D vector of Q or phi (in ascending order)
    intensity_array: ~numpy.ndarray
        2D array of intensities, shape[0] = number of Q or phi, shape[1] = number of wavelength
    error_array: ~numpy.ndarray
        2D array of errors, shape[0] = number of Q or phi, shape[1] = number of wavelength
    k_vec: ~numpy.ndarray
        calculated K vector
    k_error_vec: ~numpy.ndarray
        calculated K error vector

    Returns
    -------
    tuple
        normalized I(1D), normalized error(1D)

    """

    # Sanity check
    assert wl_vec.shape[0] == intensity_array.shape[1]  # wavelength as lambda
    assert x_vec.shape[0] == error_array.shape[0]  # points as Q or phi
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
