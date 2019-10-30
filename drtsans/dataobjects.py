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


class IQazimuthal(object):
    pass


class IQcrystal(object):
    pass
