from collections import OrderedDict

import functools
from collections import namedtuple, Mapping


class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super(MultiOrderedDict, self).__setitem__(key, value)
            # super().__setitem__(key, value) # in Python 3


def namedtuplefy(func, name=None):
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
    namedtuplefy.nt = None

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if namedtuplefy.nt is None:
            if isinstance(res, Mapping) is False:
                raise ValueError('Cannot namedtuplefy a non-dict')
            namedtuplefy.nt = \
                namedtuple(name or (func.__name__ + '_nt'), res.keys())
        return namedtuplefy.nt(**res)
    return wrapper
