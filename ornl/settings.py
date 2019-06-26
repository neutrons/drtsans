from __future__ import (absolute_import, division, print_function)

import random
import string
import inspect
from collections import OrderedDict
import functools
from collections import namedtuple, Mapping
from contextlib import contextmanager

import mantid
from mantid.simpleapi import RenameWorkspace
from mantid.api import AnalysisDataService
from mantid.kernel import ConfigService
from mantid.kernel.funcinspect import process_frame

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
    Decorator to endow a function with output_workspace optional keyword

    # Case 1: the decorated function receives an input workspace
    #         as first argument

    @optional_output_workspace
    def foo(input_ws, *args, *kwargs):
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
    def bar(input_file, *args, *kwargs):
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
        Target of decorator. Must return a workspace

    Returns
    -------
    Function
    """
    name = 'output_workspace'
    # Find out if the signature of the decorated function contains
    # optional "output_workspace" keyword
    parameters = dict(inspect.signature(func).parameters)
    output_workspace_parameter = parameters.get(name, None)
    name_in_signature = name in parameters

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        #
        # decorator does not change func if output_workspace is
        # a keyword of func and has a default value different than None
        #
        if output_workspace_parameter is not None and \
           output_workspace_parameter.default is not None:
            return func(*args, **kwargs)
        #
        # decorator does not change func if output_workspace is
        # a keyword of func and we are passing a value for it
        #
        if name_in_signature and kwargs.get(name, None) is not None:
            return func(*args, **kwargs)
        #
        # Case: output_workspace is not in func's signature but we are
        #       passing it in the list of optional arguments
        #
        if name_in_signature is False and name in kwargs:
            output_workspace = kwargs.pop(name)
            returned_workspace = func(*args, **kwargs)
            RenameWorkspace(returned_workspace,
                            OutputWorkspace=output_workspace)
            return returned_workspace
        #
        # Find out whether the first required parameter is a workspace
        #
        inspect_assignment = True
        if len(args) > 0:
            for workspace_type in workspace_types:
                if isinstance(args[0], workspace_type):
                    inspect_assignment = False
                    break
        #
        # Case we passed a workspace as first required argument
        #
        if inspect_assignment is False:
            returned_workspace = func(*args, **kwargs)
            if returned_workspace is not args[0]:  # avoid corner case
                RenameWorkspace(returned_workspace,
                                OutputWorkspace=args[0].name())
            return returned_workspace
        #
        # Inspect the name of the variable receiving the workspace
        #
        lhs = process_frame(inspect.currentframe().f_back)[1]
        output_workspace = None if len(lhs) != 1 else lhs[0]
        if output_workspace is None:
            raise RuntimeError('No output workspace name provided')
        #
        # Case: output_workspace is in func's signature but its value is None
        #
        if name_in_signature:
            kwargs[name] = output_workspace
            return func(*args, **kwargs)
        #
        # Case: output_workspace is not in func's signature and we are not
        #       passing as optional argument
        #
        returned_workspace = func(*args, **kwargs)
        RenameWorkspace(returned_workspace,
                        OutputWorkspace=output_workspace)
        return returned_workspace
    return wrapper
