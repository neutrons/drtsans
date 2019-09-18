import os
from mantid.api import mtd, MatrixWorkspace
from mantid.geometry import Instrument
from mantid.simpleapi import Load, ExtractMask
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs


def detector_name(ipt):
    r"""
    Name of the main detector array

    Parameters
    ----------
    ipt: str, Instrument, Workspace
        Input instrument name, Instrument instance, or Workspace

    Returns
    -------
    str
    """
    inst_to_det = {'EQSANS': 'detector1',
                   'EQ-SANS': 'detector1',
                   'BIOSANS': 'detector1',
                   'GPSANS': 'detector1',
                   'GenericSANS': 'detector1'}
    if isinstance(ipt, str):
        if ipt in mtd:
            instrument_name = mtd[ipt].getInstrument().getName()
        else:
            instrument_name = ipt
    elif isinstance(ipt, Instrument):
        instrument_name = ipt.getName()
    else:
        instrument_name = ipt.getInstrument().getName()  # assume workspace
    return inst_to_det[instrument_name]


def bank_detector_ids(input_workspace, masked=None):
    r"""
    Return the ID's for the detectors in detector banks (excludes monitors)

    Parameters
    ----------
    input_workspace: MatrixWorkspace
        Input workspace to find the detectors
    masked: None or Bool
        `None` yields all detector ID's; `True` yields all masked
        detector ID's; `False` yields all unmasked detector ID's

    Returns
    -------
    list
    """
    ws = mtd[str(input_workspace)]
    ids = ws.detectorInfo().detectorIDs()
    all = ids[ids >= 0].tolist()
    if masked is None:
        return all  # by convention, monitors have ID < 0
    else:
        instrument = ws.getInstrument()
        return [det_id for det_id in all if
                masked == instrument.getDetector(det_id).isMasked()]


def bank_detectors(input_workspace, masked=None):
    r"""
    Generator function to yield the detectors in the banks of the instrument.
    Excludes monitors

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace
        Input workspace to find the detectors
    masked: None or Bool
        `None` yields all detectors; `True` yields all masked detectors;
        `False` yields all unmasked detectectors
    Yields
    ------
    mantid.geometry.Detector
    """
    ws = mtd[str(input_workspace)]
    instrument = ws.getInstrument()
    for det_id in bank_detector_ids(ws, masked=masked):
        yield instrument.getDetector(det_id)
        det = instrument.getDetector(det_id)
        if masked is None or masked == det.isMasked():
            yield instrument.getDetector(det_id)


def masked_detectors(input_workspace, query_ids=None):
    r"""
    List of detector ID's that are masked

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace
        Input workspace to find the detectors
    query_ids: list
        Restrict the search to this list of detector ID's. If `None`, query
        all detectors.

    Returns
    -------
    list
    """
    mask_ws, det_ids = ExtractMask(input_workspace,
                                   OutputWorkspace=uwd())
    if query_ids is not None:
        det_ids = sorted(list(set(det_ids) & set(query_ids)))
    mask_ws.delete()
    return det_ids


def get_instrument(source):
    r"""
    Return the instrument object

    Parameters
    ----------
    source: PyObject
        Instrument object, MatrixWorkspace, , workspace name, file name,
        run number

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
        elif s in mtd:
            return get_instrument(mtd[s])
        else:
            try:
                i = int(s)
                return get_instrument(i)
            finally:
                pass

    dispatch = {Instrument: from_instrument, MatrixWorkspace: from_ws,
                int: from_integer, str: from_string}
    finder = [v for k, v in dispatch.items() if isinstance(source, k)][0]
    return finder(source)


def source_sample_distance(source, unit='mm', log_key=None, search_logs=True):
    r"""
    Report the distance (always positive!) between source and sample.

    If logs are not used or distance fails to be found in the logs, then
    calculate the distance using the instrument configuration file.

    Parameters
    ----------
    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number
    unit: str
        'mm' (millimeters), 'm' (meters)
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
        sl = SampleLogs(source)
        try:
            lk = set(log_keys).intersection(set(sl.keys())).pop()
            # uses the default unit [mm] if no unit is defined
            return float(sl.single_value(lk)) * mm2units[unit]
        except KeyError:
            pass

    # Calculate the distance using the instrument definition file
    instrument = get_instrument(source)
    sample = instrument.getSample()
    return abs(sample.getDistance(instrument.getSource())) * m2units[unit]


def sample_detector_distance(source, unit='mm', log_key=None,
                             search_logs=True):
    r"""
    Return the distance from the sample to the detector bank

    The function checks the logs for the distance, otherwise returns the
    minimum distance between the sample and the detectors of the bank

    Parameters
    ----------
    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number.
    unit: str
        'mm' (millimeters), 'm' (meters)
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
        sl = SampleLogs(source)
        try:
            lk = set(log_keys).intersection(set(sl.keys())).pop()
            # uses the default unit [mm] if no unit is defined
            return float(sl.single_value(lk)) * mm2units[unit]
        except KeyError:
            pass
    # Calculate the distance using the instrument definition file
    instrument = get_instrument(source)
    det = instrument.getComponentByName(detector_name(source))
    return det.getDistance(instrument.getSample()) * m2units[unit]


def source_detector_distance(source, unit='mm', search_logs=True):
    r"""
    Calculate distance between source and detector bank, in mili meters

    This functions is just the sum of functions `sample_source_distance`
    and `sample_detector_distance`

    Parameters
    ----------
    source: PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number
    unit: str
        'mm' (millimeters), 'm' (meters)
    search_logs: bool
        Report the value as the sum of the source to sample distance and
        sample to detector distance found in the logs

    Returns
    -------
    float

    """
    ssd = source_sample_distance(source, unit=unit, search_logs=search_logs)
    sdd = sample_detector_distance(source, unit=unit, search_logs=search_logs)
    return ssd + sdd