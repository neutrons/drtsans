from __future__ import (absolute_import, division, print_function)

import os
import numpy as np

from mantid.api import (Run, MatrixWorkspace)
from mantid.simpleapi import (Load, AddSampleLog)


class SampleLogs(object):
    r"""
    Log reader, a bit more pythonic
    """

    def __init__(self, other):
        self._run = self.find_run(other)
        self._ws = None

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
        log_unit = unit  # this could end up being None
        number_type = None  # this could end up being None

        def figure_log_type(val):
            for k, v in {list: 'Number Series', int: 'Number',
                         float: 'Number', str: 'String'}.items():
                if isinstance(val, k):
                    return v

        def figure_number_type(val):
            for k, v in {int: 'Int', float: 'Double'}.items():
                if isinstance(val, k):
                    return v
            if isinstance(val, list):
                return figure_number_type(val[0])

        log_type = figure_log_type(value)
        if log_type == 'Number':
            number_type = figure_number_type(value)

        # cast `value` to str
        if log_type == 'NumberSeries':
            log_text = ', '.join([str(x) for x in value])
        else:
            log_text = str(value)

        # Done, call Mantid algorithm
        kw = {k: v for (k,v) in zip(('LogUnit', 'NumberType'), (log_unit, number_type)) if v is not None}
        AddSampleLog(self._ws, LogName=name, LogText=log_text, LogType=log_type, **kw)

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
