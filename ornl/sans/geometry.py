from __future__ import (absolute_import, division, print_function)

import os
from mantid.api import MatrixWorkspace
from mantid.geometry import Instrument
from mantid.simpleapi import Load
from ornl.sans.samplelogs import SampleLogs


def bank_detector_ids(ws, masked=None):
    r"""
    Return the ID's for the detectors in detector banks (excludes monitors)

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace to find the detectors
    masked: None or Bool
        `None` yields all detector ID's; `True` yields all masked
        detector ID's; `False` yields all unmasked detector ID's

    Returns
    -------
    list
    """
    ids = ws.detectorInfo().detectorIDs()
    all = ids[ids >= 0].tolist()
    if masked is None:
        return all  # by convention, monitors have ID < 0
    else:
        instrument = ws.getInstrument()
        return [det_id for det_id in all if
                masked == instrument.getDetector(det_id).isMasked()]


def bank_detectors(ws, masked=None):
    r"""
    Generator function to yield the detectors in the banks of the instrument.
    Excludes monitors

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace to find the detectors
    masked: None or Bool
        `None` yields all detectors; `True` yields all masked detectors;
        `False` yields all unmasked detectectors
    Yields
    ------
    mantid.geometry.Detector
    """
    instrument = ws.getInstrument()
    for det_id in bank_detector_ids(ws, masked=masked):
        yield instrument.getDetector(det_id)
        det = instrument.getDetector(det_id)
        if masked is None or masked == det.isMasked():
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


def sample_source_distance(other, units='mm', log_key=None, search_logs=True):
    r"""
    Report the distance (always positive!) between sample and source.

    If logs are not used or distance fails to be found in the logs, then
    calculate the distance using the instrument configuration file.

    Parameters
    ----------
    other: PyObject
        Instrument object, MatrixWorkspace, file name, run number
    units: str
        'mm' (mili-meters), 'm' (meters)
    log_key: str
        Only search for the given string in the logs. Do not use default
        log keys
    search_logs: bool
        Report the value found in the logs.

    Returns
    -------
    float
        distance between source and sample, in selected units
    """
    m2units = dict(mm=1e3, m=1.0)
    mm2units = dict(mm=1.0, m=1e-3)

    # Search the logs for the distance
    if search_logs is True:
        log_keys = ('source-sample-distance', 'source_sample-distance',
                    'source_sample_distance', 'sample-source-distance',
                    'sample_source-distance', 'sample_source_distance')
        if log_key is not None:
            log_keys = [log_key]
        sl = SampleLogs(other)
        try:
            lk = set(log_keys).intersection(set(sl.keys())).pop()
            assert sl[lk].units == 'mm'  # required entry in mili meters
            return float(sl.single_value(lk)) * mm2units[units]
        except KeyError:
            pass

    # Calculate the distance using the instrument definition file
    instrument = get_instrument(other)
    sample = instrument.getSample()
    return abs(sample.getDistance(instrument.getSource())) * m2units[units]


def source_sample_distance(other, units='mm', log_key=None, search_logs=True):
    r"""Syntactic sugar of function `sample_source_distance`"""
    return sample_source_distance(other, units=units,
                                  log_key=log_key, search_logs=search_logs)


def sample_detector_distance(ws, units='mm', log_key=None, search_logs=True):
    r"""
    Return the distance from the sample to the detector bank

    The function checks the logs for the distance, otherwise returns the
    minimum distance between the sample and the detectors of the bank

    Parameters
    ----------
    ws: Workspace
        Instrument object, MatrixWorkspace, file name, run number
    units: str
        'mm' (mili-meters), 'm' (meters)
    log_key: str
        Only search for the given string in the logs. Do not use default
        log keys
    search_logs: bool
        Report the value found in the logs.

    Returns
    -------
    float
        distance between sample and detector, in selected units
    """
    m2units = dict(mm=1e3, m=1.0)
    mm2units = dict(mm=1.0, m=1e-3)

    # Search the logs for the distance
    if search_logs is True:
        log_keys = ('detector-sample-distance', 'detector_sample-distance',
                    'detector_sample_distance', 'sample-detector-distance',
                    'sample_detector-distance', 'sample_detector_distance')
        if log_key is not None:
            log_keys = [log_key]
        sl = SampleLogs(ws)
        try:
            lk = set(log_keys).intersection(set(sl.keys())).pop()
            assert sl[lk].units == 'mm'  # required entry in mili meters
            return float(sl.single_value(lk)) * mm2units[units]
        except KeyError:
            pass
    # Calculate the distance using the instrument definition file
    instrument = get_instrument(ws)
    sample = instrument.getSample()
    sdd_i = [det.getDistance(sample) for det in bank_detectors(ws)]
    return min(sdd_i) * m2units[units]


def source_detector_distance(ws, units='mm', search_logs=True):
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
    search_logs: bool
        Report the value as the sum of the source to sample distance and
        sample to detector distance found in the logs

    Returns
    -------
    float

    """
    return source_sample_distance(ws, units=units, search_logs=search_logs) +\
        sample_detector_distance(ws, units=units, search_logs=search_logs)
