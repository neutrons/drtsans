from __future__ import (absolute_import, division, print_function)

import random
import string
import inspect
from collections import OrderedDict
import functools
from collections import namedtuple, Mapping
from contextlib import contextmanager

import mantid
from mantid.api import AnalysisDataService
from mantid.kernel import ConfigService

# import mantid's workspace types exposed to python
workspace_types = [getattr(mantid.dataobjects, w_type_name) for w_type_name in
                   [s for s in dir(mantid.dataobjects) if "Workspace" in s]]


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


def optional_output_workspace(func):
    r"""
    Decorator to endow a function with output_workspace optional keyword. If
    there is a parameter called 'input_workspace' the decorator will add a
    default value of output_workspace to be input_workspace. Without an
    'input_workspace' the default value of 'output_workspace' is the result
    of :py:obj:`ornl.settings.unique_workspace_dundername`

    # Case 1: the decorated function receives an input workspace
    #         as first argument

    @optional_output_workspace
    def foo(input_workspace, *args, **kwargs):
        ....
        return some_workspace

    ws = foo(input_ws)
    print(ws.name())  # prints input_ws.name(), input_ws is overwritten.

    ws = foo(input_ws, output_workspace='new_ws')
    print(ws.name())  # prints "new_ws". output_workspace need not be part
                      # of foo's signature

    # Case 2: the decorated function does not receive an input workspace
    #         as first argument
    @optional_output_workspace
    def bar(input_file, *args, **kwargs):
        .....
        return some_workspace

    ws = bar(my_input_file)
    print(ws.name())  # prints "ws" ( introspection of  assignemnt statement)
    ws = bar(my_input_file, output_workspace='new_ws')
    print(ws.name())  # prints "new_ws".  output_workspace need not be part
                      # of bar's signature

    Parameters
    ----------
    func: Function
        Target of decorator. Must return a workspace and must accept kwargs

    Returns
    -------
    Function
    """
    name = 'output_workspace'
    # Find out if the signature of the decorated function contains
    # optional "output_workspace" keyword
    parameters = dict(inspect.signature(func).parameters)
    output_workspace_parameter = parameters.get(name, None)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        name_specified = output_workspace_parameter is not None and \
           output_workspace_parameter.default is not None

        if name_specified:
            return func(*args, **kwargs)

        # decorator does not change func if output_workspace is
        # a keyword of func and has a default value different than None
        if output_workspace_parameter is not None and \
           output_workspace_parameter.default is not None:
            return func(*args, **kwargs)

        # decorator does not change func if output_workspace is
        # a keyword of func and we are passing a value for it
        if kwargs.get(name, None) is not None:
            return func(*args, **kwargs)

        if 'input_workspace' in parameters:
            # copy from input to output
            kwargs[name] = str(parameters.get('input_workspace'))
            return func(*args, **kwargs)

        # create a random hidden workspace name
        kwargs[name] = unique_workspace_dundername()
        return func(*args, **kwargs)

    return wrapper
