from __future__ import (absolute_import, division, print_function)

import os
import numpy as np

from mantid.api import (Run, MatrixWorkspace)
from mantid.simpleapi import (mtd, Load)


class SampleLogs(object):
    r"""
    Log reader, a bit more pythonic

    Params
    ------
    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number
    """

    def __init__(self, source):
        self._ws = None
        self._run = self.find_run(source)

    def __getitem__(self, item):
        if item in self._run.keys():
            return self._run[item]

    def __getattr__(self, item):
        _run = self.__dict__['_run']
        try:
            return getattr(_run, item)
        except AttributeError:
            if item in _run.keys():
                return _run.getProperty(item)

    def insert(self, name, value, unit=None):
        r"""
        Wrapper to Mantid AddSampleLog algorithm

        The following properties of AddSampleLog are determined by
        inspection of `value`: LogText, LogType, NumberType

        Parameters
        ----------
        name: str
            log entry name
        value: str, int, double, list
            Value to insert
        unit: str
            Log unit
        """
        if not unit:
            unit = ''
        if isinstance(value, list):
            value = value[0]  # copies AddSampleLog behavior

        self._ws.mutableRun().addProperty(name, value, unit, True)

    @property
    def mantid_logs(self):
        return self._run

    @property
    def workspace(self):
        return self._ws

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
        def from_ws(ws):
            self._ws = ws
            return ws.getRun()

        def from_run(a_run):
            return a_run

        def from_integer(run_number):
            w = Load(Filename=str(run_number))
            return self.find_run(w)

        def from_string(s):
            # see if it is a file
            if os.path.isfile(s):
                w = Load(Filename=s)
                return self.find_run(w)
            # see if it is an already named data object
            elif s in mtd:
                return self.find_run(mtd[s])
            else:
                try:
                    i = int(s)
                    return self.find_run(i)
                finally:
                    pass

        dispatch = {Run: from_run, MatrixWorkspace: from_ws, int: from_integer,
                    str: from_string}
        # finder = [v for k, v in dispatch.items() if isinstance(other, k)][0]

        # If others is not None: raise exception
        if other is None:
            a = ('[DEBUG 187] dispatch: {}'.format(dispatch.items()))
            b = ('[DEBUG 187] other: {}'.format(other))
            raise NotImplementedError('{}\n{}'.format(a, b))

        finders = [v for k, v in dispatch.items() if isinstance(other, k)]
        if len(finders) == 0:
            # In case no items found
            raise RuntimeError('Input "other" of value {} is not supported to retrieve Mantid "run" object'.format(other))
        finder = finders[0]

        return finder(other)

    def find_log_with_units(self, log_key, unit=None):
        r"""
        Find a log entry in the logs, and ensure it has the right units

        Parameters
        ----------
        log_key: string
                 key of the log to find
        unit: None or string
               units string to enforce

        Returns
        -------
            log value
        """
        if log_key in self.keys():
            if unit is not None and not self[log_key].units == unit:
                error_msg = "Found %s with wrong units" % log_key
                error_msg += " [%s]" % self[log_key].units
                raise RuntimeError(error_msg)
            return float(self[log_key].value)
        raise RuntimeError("Could not find %s in logs" % log_key)
