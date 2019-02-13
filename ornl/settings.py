from __future__ import (absolute_import, division, print_function)

import random
import string
from collections import OrderedDict
import functools
from collections import namedtuple, Mapping
from contextlib import contextmanager

from mantid.api import AnalysisDataService
from mantid.kernel import ConfigService


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
                raise ValueError('Cannot namedtuplefy a non-dict')
            wrapper.nt = namedtuple(func.__name__ + '_nt', res.keys())
        return wrapper.nt(**res)
    wrapper.nt = None
    return wrapper


@contextmanager
def amend_config(new_config):
    r"""
    Context manager to safely modify Mantid Configuration Service while
    the function is executed.

    Parameters
    ----------
    new_config: dict
        (key, value) pairs to modify the configuration service
    """
    backup = dict()
    config = ConfigService.Instance()
    for key, val in new_config.items():
        backup[key] = config[key]
        config[key] = val  # config does not have an 'update' method
    try:
        yield
    finally:
        for key in new_config:
            config[key] = backup[key]


def unique_workspace_name(n=5, prefix='', suffix=''):
    r"""
    Create a random sequence of `n` lowercase characters that is guaranteed
    not to collide with the name of any existing Mantid workspace

    uws stands for Unique Workspace Name

    Parameters
    ----------
    n: int
        Size of the sequence
    prefix: str
        String to prefix the randon sequence
    suffix: str
        String to suffix the randon sequence

    Returns
    -------
    string
    """

    ws_name = ''.join(random.choice(string.ascii_lowercase) for _ in range(n))
    ws_name = '{}{}{}'.format(str(prefix), ws_name, str(suffix))
    while ws_name in AnalysisDataService.getObjectNames():
        characters = [random.choice(string.ascii_lowercase) for _ in range(n)]
        ws_name = ''.join(characters)
        ws_name = '{}{}{}'.format(str(prefix), ws_name, str(suffix))
    return ws_name


def unique_workspace_dundername():
    return unique_workspace_name(n=5, prefix='__')
