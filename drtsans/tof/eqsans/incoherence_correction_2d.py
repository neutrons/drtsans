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
        Boolean numpy array representing q_subset

    """
    _q_by_wavelength = i_of_q.intensity.reshape((qx_len*qy_len, wavelength_len))
    _mask_squeezed = np.all(np.isfinite(_q_by_wavelength), 1)
    return _mask_squeezed.repeat(wavelength_len)


def calculate_b2d(i_of_q, q_subset_mask, wavelength_len, min_incoh=False):
    """Calculates the 2D b parameters

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        Input reshaped I(Qx, Qy, wavelength)
    q_subset_mask: ~numpy.ndarray
        Boolean array defining q_subset
    wavelength_len: int
        Number of unique wavelength values
    min_incoh: bool
        Whether to select minimum b wavelength for b recalculation

    Returns
    -------
    ~numpy.ndarray
        calculated 2D b values

    """
    pass
