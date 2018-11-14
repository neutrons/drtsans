from __future__ import (absolute_import, division, print_function)

import functools
from collections import namedtuple


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
    nt = None

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        nonlocal nt
        if nt is None:
            nt = namedtuple(name or (func.__name__ + '_nt'), res.keys())
        return nt(**res)
    return wrapper


def find_return_dict(function_source):
    r"""
    Return the string representation of the return statement.

    It is assumed that the function returns a dictionary, and that the type
    of all dictionary keys is str

    Parameters
    ----------
    function_source: str
        Source code of a function

    Returns
    -------
    str
        String representation of a dictionary
    """
    line = function_source.split('return ')[-1].strip(' \n')
    if '\n' in line:
        new_line = ''
        for l in line.split('\n'):
            if '#' in l:
                new_line += l.split('#')[0].strip()  # remove comment
        line = new_line
    if '#' in line:
        line = line.split('#')[0].strip()  # remove last comment
    return line


def dict_keys(query):
    r"""
    Return the keys from the string representation of a dictionary

    Parameters
    ----------
    query: str
        String representation of a dictionary

    Returns
    -------
    list
        keys of the dictionary

    Raises
    ------
    ValueError
        input string is not a valid representation of a dictionary
    """
    if query[0:5] == 'dict(':
        key_vals = query[5:][:-1].split(',')
        return [key_val.split('=')[0].strip() for key_val in key_vals]
    elif query[0] == '{' and query[-1] == '}':
        key_vals = query[1:][:-1].split(',')
        return [key_val.split(':')[0].strip(' "\'') for key_val in key_vals]
    else:
        raise ValueError('No valid string representation of a "dict"')

'''
def namedtuplefy(func):
    r"""
    Decorator to transform the return dictionary of a function into
    a namedtuple

    Parameters
    ----------
    func: Function
        Function to be decorated

    Returns
    -------
    Function
    """
    keys = dict_keys(find_return_dict(inspect.getsource(func)))
    nt = namedtuple(func.__name__ + '_nametuplefied', keys)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return nt(**func(*args, **kwargs))
    return wrapper
'''