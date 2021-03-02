import numpy as np

from drtsans.dataobjects import IQazimuthal


def reshape_q_azimuthal(i_of_q):
    """Enforce usable IQazimuthal setup and preserve original shape

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
    index_sorted = np.lexsort((
        flat_i_of_q.wavelength,
        flat_i_of_q.qy,
        flat_i_of_q.qx
    ))
    kwargs = dict()
    if flat_i_of_q.delta_qx is not None:
        kwargs['delta_qx'] = flat_i_of_q.delta_qx[index_sorted]
    if flat_i_of_q.delta_qy is not None:
        kwargs['delta_qy'] = flat_i_of_q.delta_qy[index_sorted]
    return IQazimuthal(
        intensity=flat_i_of_q.intensity[index_sorted],
        error=flat_i_of_q.error[index_sorted],
        qx=flat_i_of_q.qx[index_sorted],
        qy=flat_i_of_q.qy[index_sorted],
        wavelength=flat_i_of_q.wavelength[index_sorted],
        **kwargs
    )


def gen_q_subset_mask(i_of_q, qx_len, qy_len, wavelength_len):
    """Generate filtering array for q_subset from intensity vector for 2d case

    Calculation of B requires usage of only qx, qy where all lambda exist,
    so this function handles creating a boolean array for filtering

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
        Boolean numpy array representing q_subset in form Qxy, wavelength

    """
    _q_by_wavelength = i_of_q.intensity.reshape((qx_len*qy_len, wavelength_len))
    _mask_squeezed = np.all(np.isfinite(_q_by_wavelength), 1)
    return _mask_squeezed


def calculate_b2d(i_of_q, q_subset_mask, qx_len, qy_len, wavelength_len, min_incoh=False):
    """Calculates the 2D b parameters

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
    min_incoh: bool
        Whether to select minimum b wavelength for b recalculation

    Returns
    -------
    tuple
        calculated 2D b values, calculated 2d b errors, reference wavelength index

    """
    _qxy = qx_len * qy_len
    # reshape to Qxy by wavelength, filter wavelengths, swap to wavelength by Qxy
    _sub_i = i_of_q.intensity.reshape((_qxy, wavelength_len))[q_subset_mask, :].transpose()
    _sub_i_e = i_of_q.error.reshape((_qxy, wavelength_len))[q_subset_mask, :].transpose()
    # initially calculate b2d using smallest wavelength as reference
    _ref = 0
    b2d, b2d_e = _b_math(_ref, _sub_i, _sub_i_e, wavelength_len)
    if min_incoh is False:
        return b2d, b2d_e, _ref
    # if min_incoh, redo calculation with minimum b wavelength as ref
    _ref = np.argmin(b2d)
    b2d, b2d_e = _b_math(_ref, _sub_i, _sub_i_e, wavelength_len)
    return b2d, b2d_e, _ref


def _b_math(ref, sub, sub_e, w_len):
    sub_len = sub.shape[0]
    # expand reference wavelength across wavelength values of subset of intensities
    ref_i = np.tile(sub[ref], w_len).reshape((w_len, sub_len))
    ref_i_e = np.tile(sub_e[ref], w_len).reshape((w_len, sub_len))
    c_val = 1/sub_len
    # do the actual math for b and b error
    b_val = -c_val * np.sum(ref_i - sub, 1)
    b_e_val = c_val * np.sqrt(np.sum(ref_i_e**2 + sub_e**2, 1))
    return b_val, b_e_val


def intensity_error(i_of_q, q_subset_mask, qx_len, qy_len, wavelength_len, ref):
    """Calculates corrected error from reshaped i_of_q and ref

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

    Returns
    -------
    ~numpy.ndarray
        intensity error vector

    """
    # collapse error and mask to 3D numpy arrays
    _i_e_pack = i_of_q.error.copy().reshape((qx_len, qy_len, wavelength_len))
    _mask_pack = q_subset_mask.reshape((qx_len, qy_len, wavelength_len))
    # filter only the reference wavelength
    _e_ref_pack = _i_e_pack[_mask_pack[:, :, ref], ref]
    # calculate summation from reference intensity errors
    _e_ref_term = np.sum(_e_ref_pack)/(_e_ref_pack.shape[0]**2)
    # step through each wavelength
    for _wave in range(wavelength_len):
        # grab slice of filter and calculate sum of filtered values
        _w_mask = _mask_pack[:, :, _wave]
        _e_w_pack = _i_e_pack[_w_mask, _wave]
        _e_w_term = np.sum(_e_w_pack)/(_e_w_pack.shape[0]**2)
        _c_val = 1 - 2/_e_w_pack.shape[0]
        # apply error correction to q_subset
        _i_e_pack[_w_mask, _wave] = _c_val*_i_e_pack[_w_mask, _wave]**2 + _e_w_term + _e_ref_term
        # apply error correction to not q_subset
        _i_e_pack[~_w_mask, _wave] = _i_e_pack[~_w_mask, _wave]**2 + _e_w_term + _e_ref_term
    # final sqrt to finish correction
    _i_e_pack = np.sqrt(_i_e_pack)
    # return flattened for consistency
    return _i_e_pack.flatten()
