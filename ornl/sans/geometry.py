from __future__ import (absolute_import, division, print_function)

import os
from mantid.api import MatrixWorkspace
from mantid.geometry import Instrument
from mantid.simpleapi import Load
from ornl.sans.samplelogs import SampleLogs


def get_instrument(other):
    r"""
    Return the instrument object

    Parameters
    ----------
    other: Instrument object, MatrixWorkspace, file name, run number

    Returns
    -------
    Mantid::Instrument
        Instrument object
    """
    def from_instrument(instrument):
        return instrument

    def from_ws(ws):
        return ws.getInstrument()

    def from_integer(run_number):
        w = Load(Filename=str(run_number))
        return get_instrument(w)

    def from_string(s):
        if os.path.isfile(s):
            w = Load(Filename=s)
            return get_instrument(w)
        else:
            try:
                i = int(s)
                return get_instrument(i)
            finally:
                pass

    dispatch = {Instrument: from_instrument, MatrixWorkspace: from_ws,
                int: from_integer, str: from_string}
    finder = [v for k, v in dispatch.items() if isinstance(other, k)][0]
    return finder(other)


def sample_source_distance(other, units='mm'):
    r"""
    Calculate the distance (always positive!) between sample and source
    using the instrument configuration file

    Parameters
    ----------
    other: Instrument object, MatrixWorkspace, file name, run number
    units: str
        'mm' (mili-meters), 'm' (meters)

    Returns
    -------
    float
        distance between source and sample, in selected units
    """
    scaling = dict(mm=1e3, m=1.0)
    instrument = get_instrument(other)
    sample = instrument.getSample()
    return abs(sample.getDistance(instrument.getSource())) * scaling[units]


def sample_detector_distance(ws, log_key='sample-detector-distance'):
    r"""
    Return the distance from the sample to the detector bank, in mili meters

    The function checks the logs for the distance, otherwise returns the
    minimum distance between the sample and any of the detector pixels.

    Parameters
    ----------
    ws: Matrixworkspace
        Workspace containing logs and a full instrument
    log_key: str
        Log entry containing the distance

    Returns
    -------
    float
    """
    sl = SampleLogs(ws)
    if log_key in sl.keys():
        assert sl[log_key].units == 'mm'
        return float(sl[log_key].value.mean())
    else:
        instrument = ws.getInstrument()
        sample = instrument.getSample()
        sdd_i = [instrument.getDetector(i).getDistance(sample)
                 for i in range(ws.getNumberHistograms())]
        return 1e03 * min(sdd_i)


def source_detector_distance(ws):
    r"""
    Calculate distance between source and detector bank, in mili meters

    This functions is just the sum of functions `sample_source_distance`
    and `sample_detector_distance`

    Parameters
    ----------
    ws Matrixworkspace
        Workspace containing logs and a full instrument

    Returns
    -------
    float

    """
    return sample_source_distance(ws) + sample_detector_distance(ws)
