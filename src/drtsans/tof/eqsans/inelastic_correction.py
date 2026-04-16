"""
This module provides the functionality to correct I(Q, wavelength) accounting for wavelength-dependent
inelastic scattering, for both 1D and 2D data.
"""

# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689
# Step 3: Correct wave-dependent incoherence intensity for I(Q, wavelength)
from drtsans.tof.eqsans.elastic_correction import (
    build_i1d_from_intensity_domain_meshgrid,
    build_i1d_one_wl_from_intensity_domain_meshgrid,
    determine_common_domain_range_mesh,
    determine_reference_wavelength_intensity_mesh,
    reshape_intensity_domain_meshgrid,
)
from drtsans.dataobjects import getDataType, DataType, save_i1d, IQazimuthal, IQmod
from collections import namedtuple
from mantid.kernel import logger
from typing import Tuple, List, Optional
import numpy as np
import os


__all__ = [
    "correct_incoherence_inelastic_all",
    "calculate_incoherence_correction_factors",
    "apply_incoherence_correction_to_unbinned_data",
    "correct_incoherence_inelastic_1d",
    "correct_incoherence_inelastic_2d",
    "CorrectedI1D",
    "CorrectedIQ2D",
    "CorrectionFactors",
    "inelastic_correction",
]

# Output of corrected 1D case
CorrectedI1D = namedtuple("CorrectedI1D", "i1d b_factor b_error")
# Output of corrected 2D case
CorrectedIQ2D = namedtuple("CorrectedIQ2D", "iq2d b_factor b_error")
# Output of correction factors calculation
CorrectionFactors = namedtuple("CorrectionFactors", "b_factor b_error wavelength")


# ====================================================================================
# 2D Incoherence Correction Functions
# ====================================================================================


def reshape_q_azimuthal(i_of_q):
    """Enforce flattened and sorted IQazimuthal setup

    2D incoherence correction as implemented operates on IQazimuthal data assuming that
    the numpy arrays are one dimensional and sorted such that a 3D array of each IQazimuthal
    attribute structured [Qx index, Qy index, wavelength index] would equal the desired array
    when flattened by numpy.ndarray.flatten

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        Input I(Qx, Qy, wavelength)

    Returns
    -------
    ~drtsans.dataobject.IQazimuthal
        flattened and sorted input

    """
    flat_i_of_q = i_of_q.ravel()
    # lexsort sorts from last argument to first (qx, then qy, then wavelength in this case)
    # this recreates the ordering of numpy.ndarray.flatten on array form [Qx, Qy, wavelength]
    index_sorted = np.lexsort((flat_i_of_q.wavelength, flat_i_of_q.qy, flat_i_of_q.qx))
    kwargs = dict()
    if flat_i_of_q.delta_qx is not None:
        kwargs["delta_qx"] = flat_i_of_q.delta_qx[index_sorted]
    if flat_i_of_q.delta_qy is not None:
        kwargs["delta_qy"] = flat_i_of_q.delta_qy[index_sorted]
    return IQazimuthal(
        intensity=flat_i_of_q.intensity[index_sorted],
        error=flat_i_of_q.error[index_sorted],
        qx=flat_i_of_q.qx[index_sorted],
        qy=flat_i_of_q.qy[index_sorted],
        wavelength=flat_i_of_q.wavelength[index_sorted],
        **kwargs,
    )


def correct_incoherence_inelastic_2d(i_of_q, b_array):
    """Correct I(Q2D) with wavelength dependent incoherence inelastic scattering

    This method implements the workflow for correcting I(Q2D) with
    wavelength-dependent incoherent inelastic scattering

    The correction of I(Q2D) uses the correction term b calculated from the 1D I(Q, lambda), since
    this gives a lower statistical error in b. The error in the corrected I(Q2D) is calculated
    assuming b is independent of I(Q2D).

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy, wavelength) with error
    b_array: ~numpy.ndarray
        2D numpy array for B[wavelength], B error[wavelength]

    Returns
    -------
    CorrectedIQ2D
        named tuple of corrected I(Qx, Qy, wavelength), b1d, b1d error

    """
    # coerce IQazimuthal data to desired shapes
    _i_of_q = reshape_q_azimuthal(i_of_q)

    # grab unique lengths
    _qx_len = np.unique(_i_of_q.qx).shape[0]
    _qy_len = np.unique(_i_of_q.qy).shape[0]

    # get b values
    b1d, b1d_error = b_array

    corrected_intensity = _i_of_q.intensity - np.tile(b1d, _qx_len * _qy_len)
    corrected_error = np.sqrt(_i_of_q.error**2 + np.tile(b1d_error**2, _qx_len * _qy_len))
    corrected_i_of_q = IQazimuthal(
        intensity=corrected_intensity,
        error=corrected_error,
        qx=_i_of_q.qx,
        qy=_i_of_q.qy,
        wavelength=_i_of_q.wavelength,
        delta_qx=_i_of_q.delta_qx,
        delta_qy=_i_of_q.delta_qy,
    )
    corrected = CorrectedIQ2D(iq2d=corrected_i_of_q, b_factor=b1d, b_error=b1d_error)
    return corrected


# ====================================================================================
# 1D Incoherence Correction Functions
# ====================================================================================


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


def calculate_incoherence_correction_factors(
    i1d_scalar_binned,
    select_minimum_incoherence,
    intensity_weighted=False,
    qmin=None,
    qmax=None,
    factor=None,
):
    """Calculate wavelength-dependent incoherence correction factors from scalar-binned I(Q, λ)

    This function calculates the b(λ) correction factors that are independent of binning mode.
    It should be called with scalar-binned I(Q, λ) data to ensure consistent correction factors
    regardless of whether the final output is scalar, wedge, or annular binning.

    Parameters
    ----------
    i1d_scalar_binned: ~drtsans.dataobjects.IQmod | ~drtsans.dataobjects.I1DAnnular
        Scalar-binned I(Q, wavelength) data used to calculate correction factors
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

    Returns
    -------
    CorrectionFactors
        Named tuple containing b_factor, b_error, and wavelength arrays
    """
    # Convert to mesh grid I(Q) and delta I(Q), or I(phi) for annular binning
    wl_vec, x_vec, i_array, error_array, _ = reshape_intensity_domain_meshgrid(i1d_scalar_binned)

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
    b_array, _, _ = calculate_b_factors(
        wl_vec,
        x_vec,
        i_array,
        error_array,
        select_minimum_incoherence,
        xmin_index,
        xmax_index,
        intensity_weighted=intensity_weighted,
    )

    return CorrectionFactors(b_factor=b_array[0], b_error=b_array[1], wavelength=wl_vec)


def apply_incoherence_correction_to_unbinned_data(
    i_of_q_2d,
    i_of_q_1d,
    correction_factors,
):
    """Apply incoherence correction factors to unbinned I(Q, λ) and I(Qx, Qy, λ) data

    This function applies pre-calculated b(λ) correction factors to unbinned data.
    The corrected unbinned data can then be binned in any mode (scalar, wedge, annular)
    and will produce consistent results.

    Parameters
    ----------
    i_of_q_2d: ~drtsans.dataobjects.IQazimuthal
        Unbinned I(Qx, Qy, wavelength) data
    i_of_q_1d: ~drtsans.dataobjects.IQmod
        Unbinned I(Q, wavelength) data
    correction_factors: CorrectionFactors
        Correction factors from calculate_incoherence_correction_factors()

    Returns
    -------
    tuple[IQazimuthal, IQmod]
        Corrected unbinned I(Qx, Qy, λ) and I(Q, λ)
    """
    from drtsans.dataobjects import IQmod

    # Helper function to create numpy-based interpolation with edge handling
    def _make_interp_fn(x, y):
        """
        Create a 1D interpolation function using numpy.interp with
        linear interpolation and edge handling equivalent to
        scipy.interpolate.interp1d(..., bounds_error=False,
        fill_value=(y[0], y[-1])).

        Parameters
        ----------
        x: array-like
            x coordinates (must be sorted)
        y: array-like
            y values corresponding to x

        Returns
        -------
        callable
            Interpolation function that accepts new x values
        """
        x = np.asarray(x)
        y = np.asarray(y)
        left = y[0]
        right = y[-1]

        def _interp(new_x):
            return np.interp(new_x, x, y, left=left, right=right)

        return _interp

    # Create interpolation functions for b_factor and b_error
    # Use linear interpolation, extrapolate with nearest values at edges
    b_factor_interp = _make_interp_fn(
        correction_factors.wavelength,
        correction_factors.b_factor,
    )
    b_error_interp = _make_interp_fn(
        correction_factors.wavelength,
        correction_factors.b_error,
    )

    # Apply correction to 2D unbinned data
    if i_of_q_2d is not None and len(i_of_q_2d.intensity) > 0:
        # Interpolate b factors to match wavelengths in unbinned data
        b_vals_2d = b_factor_interp(i_of_q_2d.wavelength)
        b_errs_2d = b_error_interp(i_of_q_2d.wavelength)

        # Apply correction: I_corrected = I - b(λ)
        corrected_intensity_2d = i_of_q_2d.intensity - b_vals_2d

        # Error propagation: σ²_corrected = σ²_I + σ²_b
        # Note: We treat I and b as independent because:
        # 1. b(λ) is calculated from scalar-binned I(Q, λ) at specific Q values within [qmin, qmax]
        # 2. The unbinned data being corrected consists of individual events at arbitrary (Qx, Qy)
        #    values that were NOT used in the b(λ) calculation
        # Therefore, the standard independent error propagation is appropriate here.
        # This differs from the correlation-aware error model in correct_intensity_error(),
        # which applies when correcting the SAME binned data used to calculate b(λ).
        corrected_error_2d = np.sqrt(i_of_q_2d.error**2 + b_errs_2d**2)

        corrected_i_of_q_2d = IQazimuthal(
            intensity=corrected_intensity_2d,
            error=corrected_error_2d,
            qx=i_of_q_2d.qx,
            qy=i_of_q_2d.qy,
            wavelength=i_of_q_2d.wavelength,
            delta_qx=i_of_q_2d.delta_qx,
            delta_qy=i_of_q_2d.delta_qy,
        )
    else:
        # Return input as-is if None or empty
        corrected_i_of_q_2d = i_of_q_2d

    # Apply correction to 1D unbinned data
    if i_of_q_1d is not None and len(i_of_q_1d.intensity) > 0:
        # Interpolate b factors to match wavelengths in unbinned data
        b_vals_1d = b_factor_interp(i_of_q_1d.wavelength)
        b_errs_1d = b_error_interp(i_of_q_1d.wavelength)

        # Apply correction: I_corrected = I - b(λ)
        corrected_intensity_1d = i_of_q_1d.intensity - b_vals_1d

        # Error propagation: σ²_corrected = σ²_I + σ²_b (independent errors, see comment above)
        corrected_error_1d = np.sqrt(i_of_q_1d.error**2 + b_errs_1d**2)

        corrected_i_of_q_1d = IQmod(
            intensity=corrected_intensity_1d,
            error=corrected_error_1d,
            mod_q=i_of_q_1d.mod_q,
            delta_mod_q=i_of_q_1d.delta_mod_q,
            wavelength=i_of_q_1d.wavelength,
        )
    else:
        # Return input as-is if None or empty
        corrected_i_of_q_1d = i_of_q_1d

    return corrected_i_of_q_2d, corrected_i_of_q_1d


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

    # If JSON parameter "selectMinIncoh" is true
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


def inelastic_correction(
    iq2d_unbinned: IQazimuthal,
    iq1d_unbinned: IQmod,
    num_x_bins: int,
    num_y_bins: int,
    num_q1d_bins: int,
    num_q1d_bins_per_decade: int,
    decade_on_center: bool,
    bin1d_type: str,
    log_binning: bool,
    user_qmin: Optional[float],
    user_qmax: Optional[float],
    annular_bin: float,
    wedges: List[Tuple[int, int]],
    symmetric_wedges: bool,
    weighted_errors: bool,
    select_min_incoherence: bool,
    intensity_weighted: bool,
    incoh_qmin: Optional[float],
    incoh_qmax: Optional[float],
    incoh_factor: Optional[float],
    output_dir: str,
    output_filename: str,
    raw_name: str,
) -> Tuple[IQazimuthal, IQmod]:
    """Apply inelastic/incoherence correction to unbinned I(Q) data.

    This function encapsulates all inelastic correction logic:
    1. Determines Q ranges from data if not provided by user
    2. Temporarily bins data in scalar mode to calculate b(λ) factors
    3. Saves b(λ) to file
    4. Applies b(λ) correction to unbinned sample data
    5. Returns corrected unbinned data

    Parameters
    ----------
    iq2d_unbinned : IQazimuthal
        Unbinned 2D I(Qx, Qy) sample data
    iq1d_unbinned : IQmod
        Unbinned 1D I(Q) sample data
    num_x_bins : int
        Number of Qx bins for temporary binning
    num_y_bins : int
        Number of Qy bins for temporary binning
    num_q1d_bins : int
        Number of Q bins for temporary 1D binning
    num_q1d_bins_per_decade : int
        Number of bins per decade for logarithmic binning
    decade_on_center : bool
        Whether decade boundaries are on bin centers
    bin1d_type : str
        Type of 1D binning ('scalar', 'wedge', 'annular')
    log_binning : bool
        Whether to use logarithmic binning
    user_qmin : float, optional
        Minimum Q value for binning (if None, determined from data)
    user_qmax : float, optional
        Maximum Q value for binning (if None, determined from data)
    annular_bin : float
        Width of annular bin in degrees
    wedges : list
        List of (angle_min, angle_max) tuples for wedges
    symmetric_wedges : bool
        Whether to add symmetric wedges
    weighted_errors : bool
        Whether to use error-weighted binning
    select_min_incoherence : bool
        Flag to determine correction B by minimum incoherence
    intensity_weighted : bool
        Whether to use intensity-weighted B factor calculation
    incoh_qmin : float, optional
        Minimum Q for incoherence correction calculation
    incoh_qmax : float, optional
        Maximum Q for incoherence correction calculation
    incoh_factor : float, optional
        Factor for automatic Q range determination
    output_dir : str
        Output directory for correction files
    output_filename : str
        Base filename for output
    raw_name : str
        Prefix for correction file names

    Returns
    -------
    tuple
        (corrected_iq2d_unbinned, corrected_iq1d_unbinned)
    """
    # Import bin_all here to avoid circular import
    from drtsans.iq import bin_all
    from drtsans.tof.eqsans.correction_api import save_b_factor

    logger.notice("Applying inelastic/incoherent correction")

    # Determine Q ranges from data if not provided by user
    qmin = user_qmin if user_qmin is not None else iq1d_unbinned.mod_q.min()
    qmax = user_qmax if user_qmax is not None else iq1d_unbinned.mod_q.max()
    qxrange = (np.min(iq2d_unbinned.qx), np.max(iq2d_unbinned.qx))
    qyrange = (np.min(iq2d_unbinned.qy), np.max(iq2d_unbinned.qy))

    # Temporarily bin data in scalar mode to calculate b(λ) factors
    # These binned results are ONLY used for factor calculation, then discarded
    iq2d_temp_binned, iq1d_temp_binned = bin_all(
        iq2d_unbinned,
        iq1d_unbinned,
        num_x_bins,
        num_y_bins,
        n1dbins=num_q1d_bins,
        n1dbins_per_decade=num_q1d_bins_per_decade,
        decade_on_center=decade_on_center,
        bin1d_type="scalar" if bin1d_type == "wedge" else bin1d_type,
        log_scale=log_binning,
        qmin=qmin,
        qmax=qmax,
        qxrange=qxrange,
        qyrange=qyrange,
        annular_angle_bin=annular_bin,
        wedges=wedges,
        symmetric_wedges=symmetric_wedges,
        error_weighted=weighted_errors,
        n_wavelength_bin=None,
    )

    if len(iq1d_temp_binned) != 1:
        raise NotImplementedError("Expected exactly one IQmod from temporary binning")

    # Calculate b(λ) correction factors from scalar-binned I(Q, λ)
    logger.notice("Calculating inelastic/incoherent correction factors from scalar-binned I(Q, lambda)")
    correction_factors = calculate_incoherence_correction_factors(
        iq1d_temp_binned[0],
        select_min_incoherence,
        intensity_weighted,
        incoh_qmin,
        incoh_qmax,
        incoh_factor,
    )

    # Save b(λ) to file
    inelastic_output_dir = os.path.join(output_dir, "info", "inelastic_incoh", output_filename)
    os.makedirs(inelastic_output_dir, exist_ok=True)

    WavelengthContainer = namedtuple("WavelengthContainer", ["wavelength"])
    wl_container = WavelengthContainer(wavelength=correction_factors.wavelength)

    save_b_factor(
        CorrectedI1D(wl_container, correction_factors.b_factor, correction_factors.b_error),
        os.path.join(inelastic_output_dir, f"{output_filename}_inelastic_b1d_{raw_name}.dat"),
    )

    # Apply b(λ) to unbinned sample data
    logger.notice("Applying inelastic/incoherent correction to unbinned data")
    iq2d_corrected, iq1d_corrected = apply_incoherence_correction_to_unbinned_data(
        iq2d_unbinned, iq1d_unbinned, correction_factors
    )

    return iq2d_corrected, iq1d_corrected
