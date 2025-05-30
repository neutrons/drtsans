# third party imports
import mantid


# standard library imports
from collections import namedtuple, OrderedDict
from collections.abc import Mapping
import functools

# import mantid's workspace types exposed to python
workspace_types = [
    getattr(mantid.dataobjects, w_type_name)
    for w_type_name in [s for s in dir(mantid.dataobjects) if "Workspace" in s]
]


class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super(MultiOrderedDict, self).__setitem__(key, value)
            # super().__setitem__(key, value) # in Python 3


def namedtuplefy(func):
    r"""
    Decorator to transform the return dictionary of a function into
    a namedtuple

    Parameters
    ----------
    func: Function
        Function to be decorated
    name: str
        Class name for the namedtuple. If None, the name of the function
        will be used
    Returns
    -------
    Function
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if wrapper.nt is None:
            if isinstance(res, Mapping) is False:
                raise ValueError("Cannot namedtuplefy a non-dict")
            wrapper.nt = namedtuple(func.__name__ + "_nt", res.keys())
        return wrapper.nt(**res)

    wrapper.nt = None
    return wrapper


def unpack_v3d(functor, index):
    """Retain only the cartesian coordinates of the V3D object returned by ```functor```

    This function reduces the memory imprint, from a V3D object to a mere 3-component list.
    Speeds up execution by avoiding crowding the heap when interating over the detectors.
    e.g.
    x = [detectorInfo().position(i) for i in range(number_detectors)]  # number_detectors V3D objects in the heap
    x = [unpackV3D(detectorInfo.position, i) for i in range(number_detectors)]  # 100 times faster

    Parameters
    ----------
    functor: function
        Callable receiving argument ```index``` and returning a V3D object.
    index: int
        DetectorInfo, ComponentInfo, or SpectrumInfo index

    Returns
    -------
    list
    """
    v3d_vector = functor(index)
    return [v3d_vector.X(), v3d_vector.Y(), v3d_vector.Z()]
