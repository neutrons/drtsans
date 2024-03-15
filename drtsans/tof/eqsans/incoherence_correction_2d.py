from collections import namedtuple

import numpy as np

from drtsans.dataobjects import IQazimuthal


__all__ = ["correct_incoherence_inelastic_2d", "CorrectedIQ2D"]


CorrectedIQ2D = namedtuple("CorrectedIQ2D", "iq2d b_factor b_error")


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


def gen_q_subset_mask(i_of_q, qx_len, qy_len, wavelength_len):
    """Generate filtering array for q_subset from intensity vector for 2d case

    Calculation of B requires usage of only qx, qy where all lambda exist,
    so this function handles creating a boolean array for filtering. Said
    array can be described by Q_valid(Qx, Qy) flattened to a 1D array of
    length qx_len*qy_len

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        Input reshaped I(Qx, Qy, wavelength)
    qx_len: int
        Number of unique Qx values
    qy_len: int
        Number of unique Qy values
    wavelength_len: int
        Number of unique wavelength values

    Returns
    -------
    ~numpy.ndarray
        1D boolean numpy array representing q_subset in form Qxy, wavelength

    """
    _q_by_wavelength = i_of_q.intensity.reshape((qx_len * qy_len, wavelength_len))
    _mask_squeezed = np.all(np.isfinite(_q_by_wavelength), 1)
    return _mask_squeezed


def intensity_error(i_of_q, q_subset_mask, qx_len, qy_len, wavelength_len, ref, b_error):
    """Calculates corrected error from reshaped i_of_q and b error vector

    Calculation of error for q_subset is described by
    (dF_{x,y}^{lambda_i})^2 = (1-2/N_q)(dI_{x,y}^{lambda_i})^2 + (dB^{lambda_i})^2
    Calculation of error for q not in q_subset is described by
    (dF_{x,y}^{lambda_i})^2 = (dI_{x,y}^{lambda_i})^2 + (dB^{lambda_i})^2
    Calculation is described in greater detail
    https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        Input reshaped I(Qx, Qy, wavelength)
    q_subset_mask: ~numpy.ndarray
        Boolean array defining q_subset
    qx_len: int
        Number of unique Qx values
    qy_len: int
        Number of unique Qy values
    wavelength_len: int
        Number of unique wavelength values
    ref: int
        Index of reference wavelength
    b_error: ~numpy.ndarray
        Error in b factor

    Returns
    -------
    ~numpy.ndarray
        intensity error vector

    """
    # collapse error and mask to 3D numpy arrays
    _i_e_pack = i_of_q.error.copy().reshape((qx_len, qy_len, wavelength_len))
    _mask_pack = q_subset_mask.repeat(wavelength_len).reshape((qx_len, qy_len, wavelength_len))
    _num_q_subset = q_subset_mask.sum()
    _c_val = 1 - 2 / _num_q_subset
    # step through each wavelength
    for _wave in range(wavelength_len):
        # no correction for reference wavelength
        if _wave == ref:
            continue
        # grab slice of filter and calculate sum of filtered values
        _w_mask = _mask_pack[:, :, _wave]
        # apply error correction to q_subset
        _i_e_pack[_w_mask, _wave] = _c_val * _i_e_pack[_w_mask, _wave] ** 2 + b_error[_wave] ** 2
        # apply error correction to not q_subset
        _i_e_pack[~_w_mask, _wave] = _i_e_pack[~_w_mask, _wave] ** 2 + b_error[_wave] ** 2
        # final sqrt to finish correction
        _i_e_pack[:, :, _wave] = np.sqrt(_i_e_pack[:, :, _wave])
    # return flattened for consistency
    return _i_e_pack.flatten()


def correct_incoherence_inelastic_2d(i_of_q, b_array, ref_wl_index):
    """Correct I(Q2D) with wavelength dependent incoherence inelastic scattering

    This method implements the workflow for correcting I(Q2D) with
    wavelength-dependent incoherent inelastic scattering

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy, wavelength) with error
    b_array: ~numpy.ndarray
        2D numpy array for B[wavelength], B error[wavelength]
    ref_wl_index: int
        Index of reference wavelength in the wavelength vector

    Returns
    -------
    CorrectedIQ2D
        named tuple of corrected I(Qx, Qy, wavelength), b2d, b2d error

    """
    # coerce IQazimuthal data to desired shapes
    _i_of_q = reshape_q_azimuthal(i_of_q)

    # grab unique lengths
    _qx_len = np.unique(_i_of_q.qx).shape[0]
    _qy_len = np.unique(_i_of_q.qy).shape[0]
    _wavelength_len = np.unique(_i_of_q.wavelength).shape[0]

    # create mask for q_subset
    q_subset = gen_q_subset_mask(_i_of_q, _qx_len, _qy_len, _wavelength_len)

    # get b values
    b2d, b2d_error = b_array

    corrected_intensity = _i_of_q.intensity - np.tile(b2d, _qx_len * _qy_len)
    corrected_error = intensity_error(_i_of_q, q_subset, _qx_len, _qy_len, _wavelength_len, ref_wl_index, b2d_error)
    corrected_i_of_q = IQazimuthal(
        intensity=corrected_intensity,
        error=corrected_error,
        qx=_i_of_q.qx,
        qy=_i_of_q.qy,
        wavelength=_i_of_q.wavelength,
        delta_qx=_i_of_q.delta_qx,
        delta_qy=_i_of_q.delta_qy,
    )
    corrected = CorrectedIQ2D(iq2d=corrected_i_of_q, b_factor=b2d, b_error=b2d_error)
    return corrected
