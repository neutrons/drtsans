from collections import namedtuple
import numpy as np

__all__ = ['IQmod', 'IQazimuthal', 'IQcrystal']


def _check_parallel(*args):
    '''This makes sure that all input arrays are parallel to each
    other. It assumes that the inputs are ndarrays.'''
    shape = args[0].shape
    for arg in args[1:]:
        if arg.shape != shape:
            raise TypeError('Shape mismatch ({} != {})'.format(shape, arg.shape))


class IQmod(namedtuple('IQmod', 'intensity error mod_q delta_mod_q wavelength')):
    def __new__(cls, intensity, error, mod_q, delta_mod_q=None, wavelength=None):
        # these conversions do nothing if the supplied information is already a numpy.ndarray
        intensity = np.array(intensity)
        error = np.array(error)
        mod_q = np.array(mod_q)

        # if intensity is 1d, then everything else will be if they are parallel
        if len(intensity.shape) != 1:
            raise TypeError('"intensity" must be a 1-dimensional array, found shape={}'.format(intensity.shape))

        # check that the manditory fields are parallel
        _check_parallel(intensity, error, mod_q)

        # work with optional fields
        if delta_mod_q is not None:
            delta_mod_q = np.array(delta_mod_q)
            _check_parallel(intensity, delta_mod_q)
        if wavelength is not None:
            wavelength = np.array(wavelength)
            _check_parallel(intensity, wavelength)

        # pass everything to namedtuple
        return super(IQmod, cls).__new__(cls, intensity, error, mod_q, delta_mod_q, wavelength)


class IQazimuthal(namedtuple('IQazimuthal', 'intensity error qx qy delta_qx delta_qy wavelength')):
    def __new__(cls, intensity, error, qx, qy, delta_qx=None, delta_qy=None, wavelength=None):
        # these conversions do nothing if the supplied information is already a numpy.ndarray
        intensity = np.array(intensity)
        error = np.array(error)
        qx = np.array(qx)
        qy = np.array(qy)

        # check that the manditory fields are parallel
        if len(intensity.shape) == 1:
            _check_parallel(intensity, error, qx, qy)
        elif len(intensity.shape) == 2:
            if len(qx.shape) == 1:
                _check_parallel(intensity, error)
                if intensity.shape[0] != qx.shape[0]:
                    raise TypeError('Incompatible dimensions intensity[{}] and qx[{}]'.format(intensity.shape,
                                                                                              qx.shape[0]))
                if intensity.shape[1] != qy.shape[0]:
                    raise TypeError('Incompatible dimensions intensity[{}] and qy[{}]'.format(intensity.shape,
                                                                                              qy.shape[0]))
            elif len(qx.shape) == 2:
                _check_parallel(intensity, error, qx, qy)
            else:
                raise TypeError('Qx can only be of dimension 1 or 2, found {}'.format(len(qx.shape)))
        else:
            raise TypeError('intensity can only be of dimension 1 or 2, found {}'.format(len(intensity.shape)))

        # work with optional fields
        if np.logical_xor(delta_qx is None, delta_qy is None):
            raise TypeError('Must specify either both or neither of delta_qx and delta_qy')
        if delta_qx is not None:
            delta_qx = np.array(delta_qx)
            delta_qy = np.array(delta_qy)
            _check_parallel(intensity, delta_qx, delta_qy)
        if wavelength is not None:
            wavelength = np.array(wavelength)
            _check_parallel(intensity, wavelength)

        # pass everything to namedtuple
        return super(IQazimuthal, cls).__new__(cls, intensity, error, qx, qy, delta_qx, delta_qy, wavelength)


class IQcrystal(namedtuple('IQazimuthal', 'intensity error qx qy qz delta_qx delta_qy delta_qz wavelength')):
    def __new__(cls, intensity, error, qx, qy, qz, delta_qx=None, delta_qy=None, delta_qz=None, wavelength=None):
        # these conversions do nothing if the supplied information is already a numpy.ndarray
        intensity = np.array(intensity)
        error = np.array(error)
        qx = np.array(qx)
        qy = np.array(qy)
        qz = np.array(qz)

        # check that the manditory fields are parallel
        if len(intensity.shape) != 1:
            raise NotImplementedError('Do not currently support dimension != 1, found {}'.format(len(intensity.shape)))
        _check_parallel(intensity, error)
        _check_parallel(intensity, qx, qy, qz)  # TODO make more generic

        # work with optional fields
        count = 0
        if delta_qx is not None:
            count += 1
        if delta_qy is not None:
            count += 1
        if delta_qz is not None:
            count += 1
        if not (count == 0 or count == 3):
            raise TypeError('Must specify either all or none of delta_qx, delta_qy, delta_qz')
        if delta_qx is not None:
            delta_qx = np.array(delta_qx)
            delta_qy = np.array(delta_qy)
            delta_qz = np.array(delta_qz)
            _check_parallel(intensity, delta_qx, delta_qy, delta_qz)
        if wavelength is not None:
            wavelength = np.array(wavelength)
            _check_parallel(intensity, wavelength)

        # pass everything to namedtuple
        return super(IQcrystal, cls).__new__(cls, intensity, error, qx, qy, qz,
                                             delta_qx, delta_qy, delta_qz, wavelength)
