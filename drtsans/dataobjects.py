from collections import namedtuple
from collections.abc import Iterable
from drtsans.settings import unique_workspace_dundername as uwd
from enum import Enum
# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
from mantid.simpleapi import mtd, CreateWorkspace
import numpy as np

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


def _nary_operation(iq_objects, operation, unpack=True):
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

    Returns
    -------
    ~drtsans.dataobjects.IQmod, ~drtsans.dataobjects.IQazimuthal, or ~drtsans.dataobjects.IQcrystal
    """
    assert len(set([type(iq_object) for iq_object in iq_objects])) == 1  # check all objects of same type
    new_components = list()
    for i in range(len(iq_objects[0])):  # iterate over the IQ object components
        i_components = [iq_object[i] for iq_object in iq_objects]  # collect the ith components of each object
        if True in [i_component is None for i_component in i_components]:  # is any of these None?
            new_components.append(None)
        elif unpack is True:
            new_components.append(operation(*i_components))
        else:
            new_components.append(operation(i_components))
    return iq_objects[0].__class__(*new_components)


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


class IQmod(namedtuple('IQmod', 'intensity error mod_q delta_mod_q wavelength')):
    '''This class holds the information for I(Q) scalar. All of the arrays must be 1-dimensional
    and parallel (same length). The ``delta_mod_q`` and ``wavelength`` fields are optional.'''
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

    def extract(self, selection):
        r"""
        Extract a subset of data points onto a new ~drtsans.dataobjects.IQmod object.

        Examples:
        - IQmod().extract(42)  # extract data point number 42
        - IQmod().extract(slice(None, None, 2))  # extract every other data point
        - IQmod().extract(IQmod().mod_q < 0.1)  # extract points with Q < 0.1

        Parameters
        ----------
        selection: int, slice, :ref:`~numpy.ndarray`
            Any selection that can be passed onto a :ref:`~numpy.ndarray`

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

    def id(self):
        return DataType.IQ_MOD

    def toWorkspace(self, name=None):
        # create a name if one isn't provided
        if name is None:
            name = uwd()

        dq = self.delta_mod_q
        if dq is None:
            dq = self.mod_q * 0.
        return CreateWorkspace(DataX=self.mod_q, DataY=self.intensity, DataE=self.error,
                               UnitX='momentumtransfer', OutputWorkspace=name, Dx=dq,
                               EnableLogging=False)


class IQazimuthal(namedtuple('IQazimuthal', 'intensity error qx qy delta_qx delta_qy wavelength')):
    '''This class holds the information for the azimuthal projection, I(Qx, Qy). The resolution terms,
    (``delta_qx``, ``delta_qy``) and ``wavelength`` fields are optional.

    All of the arrays must be 1-dimensional or 2-dimensional and matching length. For the 1-dimensional
    case, all of the arrays must be parallel (same length). For the 2-dimensional case, (``intensity``,
    ``error``, ``delta_qx``, ``delta_qy``, ``wavelength``) must all be parallel. However, for (``qx``,
    ``qy``), they must either (both) be 2-dimensional and parallel, or (both) 1-dimensional with
    ``len(qx) == intensity.shape[0]`` and ``len(qy) == intensity.shape[1]``.'''
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

    def id(self):
        return DataType.IQ_AZIMUTHAL


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
