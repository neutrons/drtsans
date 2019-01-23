from __future__ import (absolute_import, division, print_function)

import random
import string
from collections import OrderedDict
import functools
from collections import namedtuple, Mapping
from contextlib import contextmanager

from mantid.api import AnalysisDataService
from mantid.kernel import ConfigService
from mantid.simpleapi import Load


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


def unique_workspace_name(n=5):
    r"""
    Create a random sequence of `n` lowercase characters that is guaranteed
    not to collide with the name of any existing Mantid workspace

    uws stands for Unique Workspace Name

    Parameters
    ----------
    n: int
        Size of the sequence

    Returns
    -------
    string
    """

    ws_name = ''.join(random.choice(string.ascii_lowercase) for _ in range(n))
    while ws_name in AnalysisDataService.getObjectNames():
        characters = [random.choice(string.ascii_lowercase) for _ in range(n)]
        ws_name = ''.join(characters)
    return ws_name


def load_run(run_number, instrument, name=None, load_kwargs=None):
    r"""
    Load a run number for a given instrument, using the archive and
    temporarily setting Mantid's default instrument.

    This function calls Mantid's `Load` algorithm.

    Parameters
    ----------
    run_number: str, int
        Run number (e.g 98234, BIOSANS_98234)
    instrument: str
        Name of the instrument (e.g BIOSANS, GPSANS, EQSANS)
    name: str
        Name of the output workspace. If `None`, a random name is produced
    load_kwargs: dict
        Aditional optional arguments to Mantid's Load algorithm.

    Returns
    -------
    EventsWorkspace
        Workspace loaded from the run
    """
    if name is None:
        name = unique_workspace_name()
    with amend_config({'instrumentName': instrument,
                       'datasearch.searcharchive': 'on'}):
        return Load(Filename=run_number, OutputWorkspace=name)
