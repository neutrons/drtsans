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
    ~numpy.ndarray
        calculated 2D b values

    """
    _qxy = qx_len * qy_len
    # reshape to Qxy by wavelength, filter wavelengths, swap to wavelength by Qxy
    _sub_i = i_of_q.intensity.reshape((_qxy, wavelength_len))[q_subset_mask, :].transpose()
    _sub_len = _sub_i[0].shape[0]
    _c_val = -1/_sub_i[0].shape[0]
    # initially calculate b2d using smallest wavelength as reference
    # expand reference wavelength across wavelength values of subset of intensities
    _ref_i = np.tile(_sub_i[0], wavelength_len).reshape((wavelength_len, _sub_len))
    b2d = _c_val * np.sum(_ref_i - _sub_i, 1)
    if min_incoh is False:
        return b2d
    min_b = np.argmin(b2d)
    _ref_i = np.tile(_sub_i[min_b], wavelength_len).reshape((wavelength_len, _sub_len))
    b2d = _c_val * np.sum(_ref_i - _sub_i, 1)
    return b2d
