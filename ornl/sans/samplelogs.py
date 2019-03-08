from __future__ import (absolute_import, division, print_function)

import os
import numpy as np

from mantid.api import Run, MatrixWorkspace
from mantid.simpleapi import Load


class SampleLogs(object):
    r"""
    Thin wrapper around mantid.api.Run for more pythonic
    gettters and setters
    """

    def __init__(self, other):
        self.__dict__['_run'] = self.find_run(other)

    def __getitem__(self, item):
        if item in self._run.keys():
            return self._run[item]

    def __setitem__(self, key, value):
        _run = self.__dict__['_run']
        _run.addProperty(key, value, True)

    def __getattr__(self, item):
        _run = self.__dict__['_run']
        try:
            return getattr(_run, item)
        except AttributeError:
            if item in _run.keys():
                return _run.getProperty(item)

    def __setattr__(self, key, value):
        _run = self.__dict__['_run']
        if key == '_run':
            _run = value
        _run.addProperty(key, value, True)

    @property
    def mantid_logs(self):
        return self._run

    def single_value(self, log_key, operation=np.mean):
        _run = self.__dict__['_run']
        return float(operation(_run[log_key].value))

    def find_run(self, other):
        r"""
        Retrieve the Run object

        Parameters
        ----------
        other: Run, MatrixWorkspace, file name, run number

        Returns
        -------
        Run
            Reference to the run object
        """
        def from_run(a_run):
            return a_run

        def from_ws(ws):
            return ws.getRun()

        def from_integer(run_number):
            w = Load(Filename=str(run_number))
            return self.find_run(w)

        def from_string(s):
            if os.path.isfile(s):
                w = Load(Filename=s)
                return self.find_run(w)
            else:
                try:
                    i = int(s)
                    return self.find_run(i)
                finally:
                    pass

        dispatch = {Run: from_run, MatrixWorkspace: from_ws, int: from_integer,
                    str: from_string}
        finder = [v for k, v in dispatch.items() if isinstance(other, k)][0]
        return finder(other)

    def find_log_with_units(self, log_key, units=None):
        r"""
        Find a log entry in the logs, and ensure it has the right units

        Parameters
        ----------
        log_key: string
                 key of the log to find
        units: None or string
               units string to enforce

        Returns
        -------
            log value
        """
        if log_key in self.keys():
            if units is not None and not self[log_key].units == units:
                error_msg = "Found %s with wrong units" % log_key
                error_msg += " [%s]" % self[log_key].units
                raise RuntimeError(error_msg)
            return float(self[log_key].value)
        raise RuntimeError("Could not find %s in logs" % log_key)
