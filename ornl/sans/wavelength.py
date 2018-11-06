from __future__ import (absolute_import, division, print_function)

from sortedcontainers import SortedList
from functools import singledispatch, update_wrapper


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

    def __mul__(self, other):
        """
        Find the intersection band between two bands

        Parameters
        ----------
        band: Wband

        Returns
        -------
        Wband or None
            None if no common band or if common band is only a point
            as in the common band between Wband(0, 1) and Wband(1, 2)
        """
        if isinstance(other, Wband):
            a = self._min if self._min > other._min else other._min
            b = self._max if self._max < other._max else other._max
            if a >= b:
                return None
            return Wband(a, b)
        elif isinstance(other, Wbands):
            return other * self

    def __eq__(self, other):
        return self._min == other._min and self._max == other._max

    def __lt__(self, other):
        return self._min < other._min


class Wbands(object):
    r"""
    A list of non-overlapping wavelength bands. Useful for defining
    the bands transmitted by a set of disk choppers

    Parameters
    ----------
    args: Wband, or Wbands, or a list of Wband or Wbands
        A wavelength band or another Wbands object
    """

    def __init__(self, *args):
        self._bands = SortedList()
        self += args

    def __len__(self):
        return len(self._bands)

    def __eq__(self, other):
        return self._bands == other._bands

    def __getitem__(self, item):
        return self._bands[item]

    def _valid_add(self, band, index):
        r"""
        Check if new band intersect with existing bands.

        Parameters
        ----------
        band: Wband
            Candidate band to be inserted
        index: int
            index of _bands containing the immediate lower band

        Returns
        -------
        bool
            Tue if band does not intersect and can thus be added
        """
        if len(self) == 0:
            return True  # fist element goes in always
        if index < len(self):
            low_band = self._bands[index]
            if low_band * band is not None:
                return False
        if index + 1 < len(self):
            high_band = self._bands[index + 1]
            if high_band * band is not None:
                return False
        return True

    def __iadd__(self, other):
        r"""
        Insert one or more wavelength bands

        Parameters
        ----------
        other: Wband or iterable
            wavelength band to be inserted, or iterable containing
            Wband objects (e.g. another Wbands object)
        """
        if isinstance(other, Wband):
            index = self._bands.bisect_right(other)
            if self._valid_add(other, index):
                self._bands.add(other)
        else:
            try:
                for band in other:
                    self += band
            except TypeError as te:
                print('Argument is not iterable', te)
        return self

    def __mul__(self, other):
        r = Wbands()
        if isinstance(other, Wband):
            for band in self._bands:
                i = band * other
                if i is not None:
                    r += i  # add a Wband object
        elif isinstance(other, Wbands):
            for band in other._bands:
                i = self * band
                if i is not None:
                    r += i  # add from a Wbands object
        if len(r) == 0:
            return None
        return r
