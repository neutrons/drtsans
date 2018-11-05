from __future__ import (absolute_import, division, print_function)


class Wband(object):
    r"""A wavelength band, useful for defining the transmission band
    of a disk chopper

    Parameters
    ----------
    w_min: float
        Lower boundary wavelength
    w_max: float
        Upper boundary wavelength

    Raises
    -------
    ValueError
        Negative input values or lower boundary bigger than upper one
    """
    def __init__(self, w_min, w_max):
        if w_min < 0 or w_max < 0 or w_min > w_max:
            raise ValueError('Invalid wavelength band')
        self._min = w_min
        self._max = w_max

    @property
    def width(self):
        return self._max - self._min

    def intersect(self, band):
        """
        Find the common band between two bands

        Parameters
        ----------
        band: Wband

        Returns
        -------
        Wband or None
            None if no common band or if common band is only a point
            as in the common band between Wband(0, 1) and Wband(1, 2)
        """
        a = self._min if self._min > band._min else band._min
        b = self._max if self._max < band._max else band._max
        if a >= b:
            return None
        return Wband(a, b)
