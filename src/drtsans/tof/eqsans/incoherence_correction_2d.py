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
