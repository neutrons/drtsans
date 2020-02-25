from collections import namedtuple
from collections.abc import Iterable
from enum import Enum
import numpy as np
import pandas as pd

# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
from mantid.simpleapi import mtd, CreateWorkspace

from drtsans.settings import unique_workspace_dundername as uwd

__all__ = ['getDataType', 'DataType', 'IQmod', 'IQazimuthal', 'IQcrystal']


class DataType(Enum):
    WORKSPACE2D = 'Workspace2D'
    IQ_MOD = 'IQmod'
    IQ_AZIMUTHAL = 'IQazimuthal'
    IQ_CRYSTAL = 'IQcrystal'


def getDataType(obj):
    try:
        return DataType(obj.id())
    except AttributeError:
        name = str(obj)
        if name not in mtd:
            raise ValueError('Do not know how to get id from: {}'.format(obj))
        return DataType(mtd[name].id())


def _check_parallel(*args):
    '''This makes sure that all input arrays are parallel to each
    other. It assumes that the inputs are ndarrays.'''
    shape = args[0].shape
    for arg in args[1:]:
        if arg.shape != shape:
            raise TypeError('Shape mismatch ({} != {})'.format(shape, arg.shape))


def _nary_operation(iq_objects, operation, unpack=True, **kwargs):
    r"""
    Carry out an operation on the component arrays for each of the IQ objects.

    Examples:
    - _nary_operation((iq_1, iq_2), numpy.append, unpack=True)
    - _nary_operation((iq_1, iq_2), numpy.concatenate, unpack=False)

    Parameters
    ----------
    iq_objects: list
        A list of ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, or ~drtsans.dataobjects.IQcrystal
        objects.
    operation: function
        A function operating on a list of :ref:`~numpy.ndarray` objects, e.g. numpy.concatenate((array1, array2))
    unpack: bool
        If set to :py:obj:`True`, then ``operation`` receives an unpacked list of arrays. If set to :py:obj:`False`,
        then ``operation`` receives the list of arrays as a single argument.
        Examples: numpy.append(*(array1, array2)) versus numpy.concatenate((array1, array2))
    kwargs: dict
        Additional options to be passed to ``operation``.

    Returns
    -------
    ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, or ~drtsans.dataobjects.IQcrystal
    """
    reference_object = iq_objects[0]
    assert len(set([type(iq_object) for iq_object in iq_objects])) == 1  # check all objects of same type
    new_components = list()
    for i in range(len(reference_object)):  # iterate over the IQ object components
        i_components = [iq_object[i] for iq_object in iq_objects]  # collect the ith components of each object
        if True in [i_component is None for i_component in i_components]:  # is any of these None?
            new_components.append(None)
        elif unpack is True:
            new_components.append(operation(*i_components, **kwargs))
        else:
            new_components.append(operation(i_components, **kwargs))
    return reference_object.__class__(*new_components)


def _extract(iq_object, selection):
    r"""
    Extract a subset of data points onto a new IQ object.

    Examples:
    - iq_object.extract(42)  # extract data point number 42
    - iq_object.extract(slice(None, None, 2))  # extract every other data point
    - iq_object.extract(IQmod().mod_q < 0.1)  # extract points with Q < 0.1

    Parameters
    ----------
    iq_object: ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, ~drtsans.dataobjects.IQcrystal
    selection: int, slice, :ref:`~numpy.ndarray`
        Any selection that can be passed onto a :ref:`~numpy.ndarray`

    Returns
    -------
    ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, ~drtsans.dataobjects.IQcrystal
    """
    component_fragments = list()
    for component in iq_object:
        if component is None:
            component_fragments.append(None)
        else:
            fragment = component.__getitem__(selection)
            if isinstance(fragment, Iterable) is False:  # selection extracts only one data point
                fragment = [fragment, ]
            component_fragments.append(fragment)
    return iq_object.__class__(*component_fragments)


def scale_intensity(iq_object, scaling):
    r"""Rescale intensity and error for one IQ object.
    Relies on fields 'intensity' and 'error' being the first two components

    Parameters
    ----------
    iq_object: ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, ~drtsans.dataobjects.IQcrystal
    scaling: float

    Returns
    -------
    ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, ~drtsans.dataobjects.IQcrystal
    """
    intensity = scaling * iq_object.intensity
    error = scaling * iq_object.error
    return iq_object.__class__(intensity, error, *[iq_object[i] for i in range(2, len(iq_object))])


def q_azimuthal_to_q_modulo(Iq):
    """ this method converts Qazimuthal to Qmodulo using
        (1) Q = sqrt(Qx**2 + Qy**2)
        (2) sigmaQ = sqrt(sigmaQx**2 + sigmaQy**2)

    Parameters:
    ----------
    Iq: IQazimuthal object

    Returns:
    -------
    Iqmod: IQmod object
    """
    qx = Iq.qx
    qy = Iq.qy
    delta_qx = Iq.delta_qx
    delta_qy = Iq.delta_qy

    mod_q = np.sqrt(qx ** 2 + qy ** 2)
    delta_mode_q = np.sqrt(delta_qx ** 2 + delta_qy ** 2)

    q_azimuthal_to_q_modulo = namedtuple('q_azimuthal_to_q_modulo', 'mod_q, delta_mod_q')
    q_azimuthal_to_q_modulo.mod_q = mod_q
    q_azimuthal_to_q_modulo.delta_mod_q = delta_mode_q

    iqmod = IQmod(intensity=Iq.intensity,
                  error=Iq.error,
                  mod_q=mod_q,
                  delta_mod_q=delta_mode_q)

    return iqmod


class IQmod(namedtuple('IQmod', 'intensity error mod_q delta_mod_q wavelength')):
    r"""This class holds the information for I(Q) scalar. All of the arrays must be 1-dimensional
    and parallel (same length). The ``delta_mod_q`` and ``wavelength`` fields are optional."""

    @staticmethod
    def read_csv(file, sep=' '):
        r"""
        Load an intensity profile into a ~drtsans.dataobjects.IQmod object.

        Required file format:
        The first row must include the names for the file columns. The order of the columns is irrelevant and
        the names of the columns must be:
        - 'intensity' for profile intensities. This column is required.
        - 'error' for uncertainties in the profile intensities. This column is required.
        - 'mod_q' for values of Q. This column is required.
        - 'delta_mod_q' for uncertainties in the Q values. This column is optional.
        - 'wavelength' This column is optional.

        Example of file contents:
            intensity error mod_q
            1000.0 89.0 0.001
            90.0 8.0 0.01
            4.7 0.9 0.1

        Usage example:
        ```
        from drtsans.mono.gpsans import IQmod
        iq = IQmod.read_csv(file_name)
        ```

        Parameters
        ----------
        file: str
            Path to input file
        sep: str
            String of length 1. Field delimiter in the input file.

        Returns
        -------
        ~drtsans.dataobjects.IQmod
        """
        frame = pd.read_csv(file, sep=sep, dtype=np.float64)
        args = [frame[label].values for label in ['intensity', 'error', 'mod_q']]
        kwargs = {label: frame[label].values for label in ['delta_mod_q', 'wavelength']
                  if label in list(frame.columns)}
        return IQmod(*args, **kwargs)

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

    def __mul__(self, scaling):
        r"""Scale intensities and their uncertainties by a number"""
        return scale_intensity(self, scaling)

    def __rmul__(self, scaling):
        return self.__mul__(scaling)

    def __truediv__(self, divisor):
        r"""Divide intensities and their uncertainties by a number"""
        return self.__mul__(1.0 / divisor)

    def extract(self, selection):
        r"""
        Extract a subset of data points onto a new ~drtsans.dataobjects.IQmod object.

        Examples:
        - IQmod().extract(42)  # extract data point number 42
        - IQmod().extract(slice(None, None, 2))  # extract every other data point
        - IQmod().extract(IQmod().mod_q < 0.1)  # extract points with Q < 0.1

        Parameters
        ----------
        selection: int, slice, ~numpy.ndarray
            Any selection that can be passed onto a ~numpy.ndarray

        Returns
        -------
        ~drtsans.dataobjects.IQmod
        """
        return _extract(self, selection)

    def concatenate(self, other):
        r"""
        Append additional data points from another ~drtsans.dataobjects.IQmod object and return the composite as a
        new ~drtsans.dataobjects.IQmod object.

        Parameters
        ----------
        other: ~drtsans.dataobjects.IQmod
            Additional data points.

        Returns
        -------
        ~drtsans.dataobjects.IQmod
        """
        return _nary_operation((self, other), np.concatenate, unpack=False)

    def sort(self, key='mod_q'):
        r"""
        Sort the data points according to one of the components of the ~drtsans.dataobjects.IQmod object.

        Parameters
        ----------
        key: str
            Component prescribing the sorting order. Default sorting is by increasing Q value.

        Returns
        -------
        ~drtsans.dataobjects.IQmod
        """
        return _extract(self, np.argsort(getattr(self, key)))

    def id(self):
        return DataType.IQ_MOD

    def to_workspace(self, name=None):
        # create a name if one isn't provided
        if name is None:
            name = uwd()

        dq = self.delta_mod_q
        if dq is None:
            dq = self.mod_q * 0.
        return CreateWorkspace(DataX=self.mod_q, DataY=self.intensity, DataE=self.error,
                               UnitX='momentumtransfer', OutputWorkspace=name, Dx=dq,
                               EnableLogging=False)

    def to_csv(self, file, sep=' ', float_format='%.6f'):
        r"""
        Write the ~drtsans.dataobjects.IQmod object into an ASCII file.

        Parameters
        ----------
        file: str
            Path to output file
        sep: str
            String of length 1. Field delimiter for the output file.
        float_format: str
            Format string for floating point numbers.
        """
        frame = pd.DataFrame({label: value for label, value in self._asdict().items() if value is not None})
        frame.to_csv(file, index=False, sep=sep, float_format=float_format)


def load_iqmod(file, sep=' '):
    r"""
    Load an intensity profile into a ~drtsans.dataobjects.IQmod object.

    Required file format:
    The first row must include the names for the file columns. The order of the columns is irrelevant and
    the names of the columns must be:
    - 'intensity' for profile intensities. This column is required.
    - 'error' for uncertainties in the profile intensities. This column is required.
    - 'mod_q' for values of Q. This column is required.
    - 'delta_mod_q' for uncertainties in the Q values. This column is optional.
    - 'wavelength' This column is optional.

    Example of file contents:
        intensity error mod_q
        1000.0 89.0 0.001
        90.0 8.0 0.01
        4.7 0.9 0.1

    Usage example:
    ```
    from drtsans.mono.gpsans import load_iqmod
    iq = load_iqmod(file_name)
    ```

    Parameters
    ----------
    file: str
        Path to input file
    sep: str
        String of length 1. Field delimiter in the input file.

    Returns
    -------
    ~drtsans.dataobjects.IQmod
    """
    return IQmod.read_csv(file, sep=sep)


def save_iqmod(iq, file, sep=' ', float_format='%.6f'):
    r"""
    Write the ~drtsans.dataobjects.IQmod object into an ASCII file.

    Current output columns
    (Line 0: ) intensity error mod_q
    Expected
    (Line 0: ) mod_q intensity error mod_q_error

    Parameters
    ----------
    iq: ~drtsans.dataobjects.IQmod
        Profile to be saved
    file: str
        Path to output file
    sep: str
        String of length 1. Field delimiter for the output file.
    float_format: str
        Format string for floating point numbers.
    """
    iq.to_csv(file, sep=sep, float_format=float_format)


class IQazimuthal(namedtuple('IQazimuthal', 'intensity error qx qy delta_qx delta_qy wavelength')):
    '''This class holds the information for the azimuthal projection, I(Qx, Qy). The resolution terms,
    (``delta_qx``, ``delta_qy``) and ``wavelength`` fields are optional.

    All of the arrays must be 1-dimensional or 2-dimensional and matching length. For the 1-dimensional
    case, all of the arrays must be parallel (same length). For the 2-dimensional case, (``intensity``,
    ``error``, ``delta_qx``, ``delta_qy``, ``wavelength``) must all be parallel. However, for (``qx``,
    ``qy``), they must either (both) be 2-dimensional and parallel, or (both) 1-dimensional with
    ``len(qx) == intensity.shape[0]`` and ``len(qy) == intensity.shape[1]``.'''
    def __new__(cls, intensity, error, qx, qy, delta_qx=None, delta_qy=None, wavelength=None):  # noqa: C901
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

    def id(self):
        return DataType.IQ_AZIMUTHAL

    def concatenate(self, other):
        r"""
        Append additional data points from another ~drtsans.dataobjects.IQazimuthal object and
        return the composite as a new ~drtsans.dataobjects.IQazimuthal object.

        Parameters
        ----------
        other: ~drtsans.dataobjects.IQazimuthal
            Additional data points.

        Returns
        -------
        ~drtsans.dataobjects.IQazimuthal
        """
        return _nary_operation((self, other), np.concatenate, unpack=False)


class IQcrystal(namedtuple('IQazimuthal', 'intensity error qx qy qz delta_qx delta_qy delta_qz wavelength')):
    '''This class holds the information for the crystallographic projection, I(Qx, Qy, Qz). All of the
    arrays must be 1-dimensional and parallel (same length). The resolution terms, (``delta_qx``,
    ``delta_qy``, ``delta_qz``) and ``wavelength`` fields are optional.'''
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

    def id(self):
        return DataType.IQ_CRYSTAL


class _Testing:
    r"""
    Mimic the numpy.testing module by applying functions of this module to the component arrays of the IQ objects
    """

    @staticmethod
    def _nary_assertion(iq_objects, assertion_function, unpack=True, **kwargs):
        r"""
        Carry out an assertion on the component arrays for each of the IQ objects.

        Examples:
        - _nary_assertion((iq_1, iq_2), numpy.append, unpack=True)
        - _nary_assertion((iq_1, iq_2), numpy.concatenate, unpack=False)

        Parameters
        ----------
        iq_objects: list
            A list of ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, or ~drtsans.dataobjects.IQcrystal
            objects.
        assertion_function: function
            A function operating on a list of :ref:`~numpy.ndarray` objects, e.g. numpy.concatenate((array1, array2))
        unpack: bool
            If set to :py:obj:`True`, then ``assertion_function`` receives an unpacked list of arrays.
            If set to :py:obj:`False`, then ``assertion_function`` receives the list of arrays as a single argument.
            Examples: numpy.append(*(array1, array2)) versus numpy.concatenate((array1, array2))
        kwargs: dict
            Additional options to be passed

        Raises
        ------
        AssertionError
        """
        reference_object = iq_objects[0]  # pick the first of the list as reference object
        assert len(set([type(iq_object) for iq_object in iq_objects])) == 1  # check all objects of same type
        for i in range(len(reference_object)):  # iterate over the IQ object components
            component_name = reference_object._fields[i]
            i_components = [iq_object[i] for iq_object in iq_objects]  # collect the ith components of each object
            if True in [i_component is None for i_component in i_components]:  # is any of these None?
                if set(i_components) == set([None]):
                    continue  # all arrays are actually None, so they are identical
                else:
                    raise AssertionError(f'field {component_name} is None for some of the iQ objects')
            elif unpack is True:
                assertion_function(*i_components, **kwargs)
            else:
                assertion_function(i_components, **kwargs)

    @staticmethod
    def assert_allclose(actual, desired, **kwargs):
        r"""Apply :ref:`~numpy.testing.assert_allclose on each component array"""
        _Testing._nary_assertion((actual, desired), assertion_function=np.testing.assert_allclose, unpack=True,
                                 **kwargs)


testing = _Testing()  # use it as if it were a module, e.g. testing.assert_allclose(iq_1, iq_2, atol=1.e-6)
