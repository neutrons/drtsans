from __future__ import (absolute_import, division, print_function)

from mantid.api import Run, MatrixWorkspace
from mantid.simpleapi import Load

import os


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

    def __getattr__(self, item):
        _run = self.__dict__['_run']
        if item in _run.keys():
            return _run.getProperty(item)

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
