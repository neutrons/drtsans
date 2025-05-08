"""
This module provides the functionality to correct I(Q, wavelength) accounting for wavelength-dependent incoherent
inelastic scattering, for both 1D and 2D data (despite its name).
"""

# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
# Step 3: Correct wave-dependent incoherence intensity for I(Q, wavelength)
from drtsans.tof.eqsans.elastic_reference_normalization import (
    build_i1d_from_intensity_domain_meshgrid,
    build_i1d_one_wl_from_intensity_domain_meshgrid,
    determine_common_domain_range_mesh,
    determine_reference_wavelength_intensity_mesh,
    reshape_intensity_domain_meshgrid,
)
from drtsans.tof.eqsans.incoherence_correction_2d import correct_incoherence_inelastic_2d
from drtsans.dataobjects import getDataType, DataType, save_i1d
from collections import namedtuple
from mantid.kernel import logger
import numpy as np
import os


__all__ = ["correct_incoherence_inelastic_all", "correct_incoherence_inelastic_1d", "CorrectedI1D"]

# Output of corrected 1D case
CorrectedI1D = namedtuple("CorrectedI1D", "i1d b_factor b_error")


def correct_incoherence_inelastic_all(
    i_of_q_2d,
    i1d,
    select_minimum_incoherence,
    intensity_weighted=False,
    qmin=None,
    qmax=None,
    factor=None,
    output_wavelength_dependent_profile=False,
    output_dir=None,
):
    """Correct I(Q2D) and I(1D) accounting wavelength dependant incoherent inelastic scattering

    This is the envelope method for the complete workflow to correct I(1D) and I(Q2D) accounting
    wavelength-dependent incoherent inelastic scattering

    Parameters
    ----------
    i_of_q_2d: ~drtsans.dataobjects.IQazimuthal or None
        I(Q2D, wavelength) and error
    i1d: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        I(1D, wavelength) and error
    select_minimum_incoherence: bool
        flag to determine correction B by minimum incoherence
    intensity_weighted: bool
        if set to true, the B factor is calculated using weighted function by intensity
    qmin: float
        manually set the qmin used for incoherent calculation
    qmax: float
        manually set the qmax used for incoherent calculation
    factor: float
        automatically determine the qmin qmax by checking the intensity profile
    output_wavelength_dependent_profile: bool
        If True then output intensity profile for each wavelength before and after b correction
    output_dir: str
        output directory for intensity profiles

    Returns
    -------
    tuple[CorrectedIQ2D, CorrectedI1D]

    """

    # Convert to mesh grid I(Q) and delta I(Q), or I(phi) for annular binning
    # x is either Q or phi
    wl_vec, x_vec, i_array, error_array, delta_x_array = reshape_intensity_domain_meshgrid(i1d)

    if qmin is not None and qmax is not None:
        xmin_index, xmax_index = np.searchsorted(x_vec, [qmin, qmax])
        xmax_index = min(xmax_index, len(x_vec) - 1)
    else:
        # determine x min and x max that exists in all I(Q, wl) or I(phi, wl)
        xmin_index, xmax_index = determine_common_domain_range_mesh(x_vec, i_array)

    if factor is not None:
        logger.notice(f"Using automated (xmin, xmax) finder with factor={factor}")
        xmin_index, xmax_index = tune_xmin(xmin_index, xmax_index, i_array, factor=factor)

    logger.notice(
        f"Incoherent correction using xmin={x_vec[xmin_index]} xmax={x_vec[xmax_index]} "
        f"with xmin_index={xmin_index}, xmax_index={xmax_index}"
    )

    # calculate B factors and errors
    b_array, ref_wl_ie, ref_wl_index = calculate_b_factors(
        wl_vec,
        x_vec,
        i_array,
        error_array,
        select_minimum_incoherence,
        xmin_index,
        xmax_index,
        intensity_weighted=intensity_weighted,
    )

    # Correct 1D
    i1d_type = getDataType(i1d)
    corrected_1d = correct_incoherence_inelastic_1d(
        wl_vec,
        x_vec,
        i_array,
        error_array,
        delta_x_array,
        b_array,
        xmin_index,
        xmax_index,
        ref_wl_ie,
        output_wavelength_dependent_profile,
        output_dir,
        i1d_type,
    )

    # Correct 2D
    if i_of_q_2d:
        corrected_2d = correct_incoherence_inelastic_2d(i_of_q_2d, b_array)
    else:
        corrected_2d = None

    return corrected_2d, corrected_1d


def correct_incoherence_inelastic_1d(
    wl_vec,
    x_vec,
    i_array,
    error_array,
    delta_x_array,
    b_array,
    xmin_index,
    xmax_index,
    ref_wl_ie,
    output_wavelength_dependent_profile=False,
    output_dir=None,
    i1d_type=DataType.IQ_MOD,
):
    """
    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        1D vector of wavelength (in ascending order)
    x_vec: ~numpy.ndarray
        1D vector of Q or phi (in ascending order)
    i_array: ~numpy.ndarray
        2D array of intensities, shape[0] = number of Q or phi, shape[1] = number of wavelength
    error_array: ~numpy.ndarray
        2D array if intensity errors
    delta_x_array: ~numpy.ndarray
        2D array of errors, shape[0] = number of Q or phi, shape[1] = number of wavelength
    b_array: ~numpy.ndarray
        incoherence inelastic correction factor B, row 0: B factor, row 1: delta B
    xmin_index: int
        index of minimum common Q or phi in x vector
    xmax_index: int
        index of maximum common Q or phi in x vector
    ref_wl_ie: ReferenceWavelength
        the reference wavelength data used to calculate B
    output_wavelength_dependent_profile: bool
        If True then output Iq for each wavelength before and after b correction
    output_dir: str
        output directory for Iq profiles
    i1d_type: DataType
        the data type of the 1D intensity profile to output

    Returns
    -------
    CorrectedI1D
        I(Q1D) or I(phi) with inelastic incoherent correction applied
    """
    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_before_b_correction.dat")
            i1d_wl = build_i1d_one_wl_from_intensity_domain_meshgrid(
                x_vec, i_array, error_array, delta_x_array, tmpwlii, i1d_type
            )
            save_i1d(i1d_wl, tmpfn)

    # correct intensities and errors
    corrected_intensities, corrected_errors = correct_intensity_error(
        wl_vec, x_vec, i_array, error_array, b_array, xmin_index, xmax_index, ref_wl_ie
    )

    # construct the output and return
    corrected_i1d = build_i1d_from_intensity_domain_meshgrid(
        wl_vec, x_vec, corrected_intensities, corrected_errors, delta_x_array, i1d_type
    )

    if output_wavelength_dependent_profile and output_dir:
        for tmpwlii, wl in enumerate(wl_vec):
            tmpfn = os.path.join(output_dir, f"IQ_{wl:.3f}_after_b_correction.dat")
            i1d_wl = build_i1d_one_wl_from_intensity_domain_meshgrid(
                x_vec, corrected_intensities, corrected_errors, delta_x_array, tmpwlii, i1d_type
            )
            save_i1d(i1d_wl, tmpfn)

    corrected = {
        "i1d": corrected_i1d,
        "b_factor": b_array[0],
        "b_error": b_array[1],
    }

    return CorrectedI1D(**corrected)


def tune_xmin(xmin_idx, xmax_idx, i_arr, factor):
    """get I(xmin_index) and I(xmax_index) from the default xmin_index and
    xmax_index if I(xmin_index) > 10 * I(xmax_index), then xmin_index
    = xmin_index + 1 repeat comparison until I(xmin_index) <= factor*
    I(xmax_index)
    """
    assert xmin_idx < xmax_idx

    imin = np.nansum(i_arr[xmin_idx])
    imax = np.nansum(i_arr[xmax_idx])
    if imin > factor * imax:
        return tune_xmin(xmin_idx + 1, xmax_idx, i_arr, factor)

    return (xmin_idx, xmax_idx)


def calculate_b_factors(
    wl_vec,
    x_vec,
    intensity_array,
    error_array,
    select_min_incoherence,
    xmin_index,
    xmax_index,
    intensity_weighted=False,
):
    """Determine reference wavelength and then calculate B factor, B error factor.

    With option select_min_incoherence, reference wavelength will be reselected according to first
    round of calculation of B

    Parameters
    ----------
    wl_vec: ~numpy.ndarray
        wavelength vector
    x_vec: ~numpy.ndarray
        Q or phi vector
    intensity_array: ~numpy.ndarray
        intenisty 2D array
    error_array: ~numpy.ndarray
        intensity error 2D array
    select_min_incoherence: bool
        flag to apply select minimum incoherence algorithm
    xmin_index: int
        index of minimum common Q or phi
    xmax_index: int
        index of maximum common Q or phi (included)
    intensity_weighted: bool
        if set to true, the B factor is calculated using weighted function by intensity
    Returns
    -------
    tuple[~numpy.ndarray, ReferenceWavelength, int]
        - row 0: B factor, row 1: delta B
        - the reference wavelength used to calculate B
        - the index of the reference wavelength in the wavelength vector
    """
    # Sanity check
    assert intensity_array.shape == error_array.shape
    assert wl_vec.shape[0] == intensity_array.shape[1]
    assert x_vec.shape[0] == error_array.shape[0]

    # determine the reference wavelength: minimum wavelength bin
    ref_wl_ie = determine_reference_wavelength_intensity_mesh(
        wl_vec, x_vec, intensity_array, error_array, xmin_index, xmax_index, 0
    )

    # calculate b(lambda_i) and delta b(lambda_i) if it is final
    b_array = calculate_b_error_b(
        wl_vec,
        intensity_array,
        error_array,
        xmin_index,
        xmax_index,
        ref_wl_ie,
        calculate_b_error=not select_min_incoherence,
        intensity_weighted=intensity_weighted,
    )

    # If JSON parameter “selectMinIncoh” is true
    if select_min_incoherence and np.argmin(b_array[0]) > 0:
        # reselect reference wavelength to the wavelength bin with minimum b
        ref_wl_index = np.argmin(b_array[0])
        # (re)determine the reference wavelengths' intensities and errors
        ref_wl_ie = determine_reference_wavelength_intensity_mesh(
            wl_vec,
            x_vec,
            intensity_array,
            error_array,
            xmin_index,
            xmax_index,
            ref_wl_index,
        )
        # (re)calculate b array
        b_array = calculate_b_error_b(
            wl_vec,
            intensity_array,
            error_array,
            xmin_index,
            xmax_index,
            ref_wl_ie,
            calculate_b_error=True,
            intensity_weighted=intensity_weighted,
        )
        """
        We would rather see where the negative values happen if it happens.
        In addition, if we later decide to discard b(lambda) from last wavelength
        bins for getting min(b[lambda]), b[lambda_max] or b[lambda_min] may have negative values.
        """
        # verify
        """
        assert (
            b_array[np.isfinite(b_array)].min() >= -1e-20
        ), f"B array has negative values: {b_array}"
        """
    else:
        ref_wl_index = 0
    return b_array, ref_wl_ie, ref_wl_index


def calculate_b_error_b(
    wl_vec,
    intensity_array,
    error_array,
    xmin_index,
    xmax_index,
    ref_wavelengths,
    calculate_b_error,
    intensity_weighted=False,
):
    """Calculate B factor and its error

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
    xmin_index: int
        index of common Q or phi min (included)
    xmax_index: int
        index of common Q or phi max (included)
    ref_wavelengths: ReferenceWavelengths
        instance of ReferenceWavelengths containing intensities and errors
    calculate_b_error: bool
        flag to calculate B factor's error
    intensity_weighted: bool
        if set to true, the B factor is calculated using weighted function by intensity

    Returns
    -------
    ~numpy.ndarray
        row 0: B factor, row 1: delta B

    """
    # Declare B factor array
    b_factor_array = np.zeros(shape=(2, len(wl_vec)), dtype="float")

    print(f"Using intensity weighted B factor calculation: {intensity_weighted}")

    from uncertainties import unumpy

    if intensity_weighted:
        # if use reference intensity to normalize the function then b value should be
        # adjusted according to the I(Q) profile on each Q

        # Equation from https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/896
        # Calculate B factors
        # b[wl] =
        # - 1/N sum_{q_k=q_min}^{q_max}(1/RefI(q_k)) * sum_{q_k=q_min}^{q_max} ([RefI(q_k) - I(q_k, wl)]/RefI(q_k))
        num_x = xmax_index + 1 - xmin_index
        # operation into a (num_x, num_wl) 2D array

        ref_intensity_vec = unumpy.uarray(
            ref_wavelengths.intensity_vec[xmin_index : xmax_index + 1].reshape((num_x, 1)),
            ref_wavelengths.error_vec[xmin_index : xmax_index + 1].reshape((num_x, 1)),
        )
        intensity_vec = unumpy.uarray(
            intensity_array[xmin_index : xmax_index + 1, :], error_array[xmin_index : xmax_index + 1, :]
        )

        b_vec = (ref_intensity_vec - intensity_vec) / ref_intensity_vec

        b_array = -1.0 * np.sum(b_vec, axis=0) * 1 / np.sum(1 / ref_intensity_vec)

        b_factor_array[0] = unumpy.nominal_values(b_array)

        # Calculate B error (delta B) as an option
        if calculate_b_error:
            # delta b(wl)^2 = 1/N^2 sum_{q_k=qmin}^{qmax} [(delta I(q_k, ref_wl))^2 + (delta I(q_k, wl))^2]
            # operation into a (num_x, num_wl) 2D array
            b_factor_array[1] = unumpy.std_devs(b_array)

    else:
        # Calculate B factors
        # b[wl] = - 1/N sum_{q_k=q_min}^{q_max} [RefI(q_k) - I(q_k, wl)]
        num_x = xmax_index + 1 - xmin_index
        # operation into a (num_x, num_wl) 2D array
        b_vec = (
            ref_wavelengths.intensity_vec[xmin_index : xmax_index + 1].reshape((num_x, 1))
            - intensity_array[xmin_index : xmax_index + 1, :]
        )
        b_factor_array[0] = -1.0 / num_x * np.sum(b_vec, axis=0)

        # Calculate B error (delta B) as an option
        if calculate_b_error:
            # delta b(wl)^2 = 1/N^2 sum_{q_k=qmin}^{qmax} [(delta I(q_k, ref_wl))^2 + (delta I(q_k, wl))^2]
            # operation into a (num_x, num_wl) 2D array
            b2_vec = (ref_wavelengths.error_vec[xmin_index : xmax_index + 1].reshape((num_x, 1))) ** 2 + (
                error_array[xmin_index : xmax_index + 1, :]
            ) ** 2
            b_factor_array[1] = 1.0 / num_x * np.sqrt(b2_vec.sum(axis=0))

    return b_factor_array


def correct_intensity_error(
    wavelength_vec,
    x_vec,
    intensity_array,
    error_array,
    b_array2d,
    xmin_index,
    xmax_index,
    ref_wl_ie,
):
    """Correct intensity and error

    Parameters
    ----------
    wavelength_vec: ~numpy.ndarray
        1D vector of wavelength (in ascending order)
    x_vec: ~numpy.ndarray
        1D vector of Q or phi (in ascending order)
    intensity_array: ~numpy.ndarray
        2D array of intensities, shape[0] = number of Q, shape[1] = number of wavelength
    error_array: ~numpy.ndarray
        2D array of errors, shape[0] = number of Q, shape[1] = number of wavelength
    xmin_index: int
        index of minimum common Q or phi (included)
    xmax_index: int
        index of maximum common Q or phi (included)
    ref_wl_ie: ReferenceWavelengths
        instance of ReferenceWavelengths containing intensities and errors
    b_array2d: ~numpy.ndarray
        2D numpy array for B[wavelength], B error[wavelength]

    Returns
    -------
    tuple
        2D arrays for (1) corrected intensities (2) corrected intensity errors

    """

    # Sanity checks
    assert intensity_array.shape == error_array.shape
    assert wavelength_vec.shape[0] == intensity_array.shape[1]
    assert x_vec.shape[0] == error_array.shape[0]
    assert len(b_array2d.shape) == 2 and b_array2d.shape[0] == 2, (
        f"Expected input B and B error but not of shape {b_array2d.shape}"
    )
    assert b_array2d.shape[1] == wavelength_vec.shape[0]

    # Init data structure
    num_common_q = xmax_index - xmin_index + 1
    corrected_intensities = np.zeros_like(intensity_array)
    corrected_errors = np.zeros_like(error_array)

    # Loop around wavelength
    for i_wl in range(wavelength_vec.shape[0]):
        # Correct intensity: I'(q, wl_i) = I(q, wl_i) - b(wl_i)
        corrected_intensities[:, i_wl] = intensity_array[:, i_wl] - b_array2d[0][i_wl]

        # Correct intensity error
        # outside q_min and q_max
        # term1[q_j] = (error[q_j, wl]^2
        term1 = error_array[:, i_wl] ** 2

        # term2 = 1/N^2 sum_{q_k=q_min:qmax}[RefError(q_k)^2 + error(q_k, wl)^2]
        # term2 is a single value
        term2 = (
            1.0
            / num_common_q**2
            * np.sum(
                ref_wl_ie.error_vec[xmin_index : xmax_index + 1] ** 2
                + error_array[xmin_index : xmax_index + 1, i_wl] ** 2
            )
        )

        # make correct term1 for those inside qmin and qmax
        # term1[q_j] = (error[q_j, wl]^2 * (1 - 2/N)
        term1[xmin_index : xmax_index + 1] *= 1 - 2 / num_common_q

        # sum
        term1 += term2

        # assign to output
        corrected_errors[:, i_wl] = term1

    return corrected_intensities, np.sqrt(corrected_errors)
