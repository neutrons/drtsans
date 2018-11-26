from __future__ import (absolute_import, division, print_function)

import os
from mantid.api import MatrixWorkspace
from mantid.geometry import Instrument
from mantid.simpleapi import Load
from ornl.sans.samplelogs import SampleLogs


def bank_ids(ws):
    r"""
    Return the ID's for the detectors in detector banks (excludes monitors)

    Parameters
    ----------
    ws: MatrixWorkspace

    Returns
    -------
    list
    """
    ids = ws.detectorInfo().detectorIDs()
    return ids[ids >= 0].tolist()  # by convention, monitors have ID < 0


def bank_detectors(ws):
    r"""
    Generator function to yield the detectors in the banks of the instrument.
    Excludes monitors

    Parameters
    ----------
    ws: MatrixWorkspace

    Yields
    ------
    mantid.geometry.Detector
    """
    instrument = ws.getInstrument()
    for det_id in bank_ids(ws):
        yield instrument.getDetector(det_id)


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


def source_sample_distance(other, units='mm'):
    r"""Syntactic sugar of function `sample_source_distance`"""
    return sample_source_distance(other, units=units)


def sample_detector_distance(ws, log_key='sample-detector-distance',
                             units='mm'):
    r"""
    Return the distance from the sample to the detector bank

    The function checks the logs for the distance, otherwise returns the
    minimum distance between the sample and the detectors of the bank

    Parameters
    ----------
    ws: Matrixworkspace
        Workspace containing logs and a full instrument
    log_key: str
        Log entry containing the distance
    units: str
        'mm' (mili-meters), 'm' (meters)

    Returns
    -------
    float
    """
    m2units = dict(mm=1e3, m=1.0)
    mm2units = dict(mm=1.0, m=1e-3)
    sl = SampleLogs(ws)
    if log_key in sl.keys():
        assert sl[log_key].units == 'mm'  # required entry in mili meters
        return float(sl[log_key].value.mean()) * mm2units[units]
    else:
        instrument = ws.getInstrument()
        sample = instrument.getSample()
        sdd_i = [det.getDistance(sample) for det in bank_detectors(ws)]
        return min(sdd_i) * m2units[units]


def source_detector_distance(ws, units='mm'):
    r"""
    Calculate distance between source and detector bank, in mili meters

    This functions is just the sum of functions `sample_source_distance`
    and `sample_detector_distance`

    Parameters
    ----------
    ws: Matrixworkspace
        Workspace containing logs and a full instrument
    units: str
        'mm' (mili-meters), 'm' (meters)

    Returns
    -------
    float

    """
    return source_sample_distance(ws, units=units) +\
        sample_detector_distance(ws, units=units)
