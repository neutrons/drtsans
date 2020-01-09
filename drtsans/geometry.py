import os

from mantid.api import MatrixWorkspace
from mantid.geometry import Instrument
from mantid.kernel import logger
from mantid.simpleapi import mtd, Load
import numpy as np

from drtsans.samplelogs import SampleLogs
from drtsans.instruments import InstrumentEnumName, instrument_enum_name
from collections import defaultdict

__all__ = ['beam_radius', 'sample_aperture_diameter', 'source_aperture_diameter']


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
    inst_to_det = defaultdict(lambda: 'detector1')
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


def main_detector_panel(source):
    r"""
    Return the main detector panel of the instrument

    Parameters
    ----------
    source: PyObject
        Instrument object, ~mantid.api.MatrixWorkspace,  ~mantid.api.IEventsWorkspace, workspace name, file path,
        run number

    Returns
    -------
    ~mantid.geometry.CompAssembly
    """
    return get_instrument(source).getComponentByName(detector_name(source))


def bank_workspace_index_range(input_workspace, component=''):
    '''
    Returns the range of workspace indices to for the named component. If no component is
    specified it is the range for the whole instrument.

    Assumptions: 1. There is one detector per spectrum 2. The detector ids are offset from
    the workspace indices 3. The lowest detector id is the first one encountered when looping
    through the spectra

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace
        Input workspace to find the detectors
    component: str
        Name of the component to get detector ids from

    Returns
    -------
    tuple
        (workspace_index_min, workspace_index_max)
    '''
    detector_ids = bank_detector_ids(input_workspace, component, None)
    detector_id_first = detector_ids.min()

    input_workspace = mtd[str(input_workspace)]
    first = None
    for i in range(input_workspace.getNumberHistograms()):
        ids = input_workspace.getSpectrum(i).getDetectorIDs()
        if len(ids) > 1:
            raise RuntimeError('do not know how to work with more than one '
                               'detector per spectrum ({})'.format(ids))
        if ids[0] == detector_id_first:
            first = i
            break
    if first is None:
        raise RuntimeError('something meaningful goes here')
    else:
        return (first, first + detector_ids.size)


def bank_detector_ids(input_workspace, component='', masked=None):
    r"""
    Return the ID's for the detectors in detector banks (excludes monitors)

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace
        Input workspace to find the detectors
    component: str
        Name of the component to get detector ids from
    masked: None or bool
        py:obj:`None` yields all detector ID's; ``True`` yields all masked
        detector ID's; ``False`` yields all unmasked detector ID's

    Returns
    -------
    ~numpy.ndarray
    """
    # the object in mantid that knows everything
    detectorInfo = mtd[str(input_workspace)].detectorInfo()
    # the full list of detector ids. The indices are parallel to this array
    ids = detectorInfo.detectorIDs()

    # which detector indices to use
    indices_to_use = np.ndarray((ids.size), dtype=bool)  # everything starts as False
    indices_to_use.fill(True)

    # sub-select the component wanted
    if component:
        componentInfo = mtd[str(input_workspace)].componentInfo()
        componentIndex = componentInfo.indexOfAny(component)
        detectorIndices = componentInfo.detectorsInSubtree(componentIndex)

        # set the the indices to only use the detectors we are interested in
        indices_to_use.fill(False)
        indices_to_use[detectorIndices] = True
    else:
        # don't use monitors
        for i in range(detectorInfo.size()):
            if detectorInfo.isMonitor(i):
                indices_to_use[i] = False

    if masked is None:
        pass
    elif detectorInfo.hasMaskedDetectors():
        for i in range(detectorInfo.size()):
            # use based on masked/not masked
            indices_to_use[i] = detectorInfo.isMasked(i) == masked
    else:
        # if there aren't masked detectors, but they are asked for, return nothing
        if masked:
            indices_to_use.fill(False)

    return ids[indices_to_use]


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
                    'sample_source-distance', 'sample_source_distance',
                    'sample_aperture_sample_distance')
        if log_key is not None:
            log_keys = (log_key)
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
            log_keys = (log_key)
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


def sample_aperture_diameter(input_workspace, unit='mm'):
    r"""
    Find and return the sample aperture diameter from the logs.

    Log keys searched are 'sample_aperture_diameter' and additional log entries for specific instruments. It is
    assumed that the units of the logged value is mm

    Parameters
    ----------
    input_workspace: :py:obj:`~mantid.api.MatrixWorkspace`
        Input workspace from which to find the aperture
    unit: str
        return aperture in requested length unit, either 'm' or 'mm'

    Returns
    -------
    float
    """
    # Additional log keys aiding in calculating the sample aperture diameter
    additional_log_keys = {InstrumentEnumName.EQSANS: ['beamslit4'],
                           InstrumentEnumName.GPSANS: [],
                           InstrumentEnumName.BIOSANS: []}
    log_keys = ['sample_aperture_diameter'] + additional_log_keys[instrument_enum_name(input_workspace)]

    sample_logs = SampleLogs(input_workspace)
    diameter = None
    for log_key in log_keys:
        if log_key in sample_logs.keys():
            diameter = sample_logs.single_value(log_key)
            break

    if diameter is None:
        raise RuntimeError('Unable to retrieve the sample aperture diameter from the logs')

    # The diameter was found using the additional logs. Insert a log for the diameter under key
    # "sample_aperture_diameter"
    if 'sample_aperture_diameter' not in sample_logs.keys():
        sample_logs.insert('sample_aperture_diameter', diameter, unit='mm')

    return diameter if unit == 'mm' else diameter / 1.e3


def source_aperture_diameter(input_workspace, unit='mm'):
    r"""
    Find and return the sample aperture diameter from the logs, or compute this quantity from other log entries.

    Log key searched is 'source_aperture_diameter'. It is assumed that the units of the logged value is mm

    Parameters
    ----------
    input_workspace: :py:obj:`~mantid.api.MatrixWorkspace`
        Input workspace from which to find the aperture
    unit: str
        return aperture in requested length unit, either 'm' or 'mm'

    Returns
    -------
    float
    """
    try:
        diameter = SampleLogs(input_workspace).single_value('source_aperture_diameter')
    except RuntimeError:
        diameter = SampleLogs(input_workspace).single_value('source_aperture_radius')
    return diameter if unit == 'mm' else diameter / 1.e3


def beam_radius(input_workspace, unit='mm'):
    """
    Calculate the beam radius impinging on the detector

    R_beam = R_sampleAp + SDD * (R_sampleAp + R_sourceAp) / SSD, where
    R_sampleAp: radius of the sample aperture,
    SDD: distance between the sample and the detector,
    R_sourceAp: radius of the source aperture,
    SSD: distance between the source and the sample.

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace, str
        Input workspace
    unit: str
        Units of the output beam radius. Either 'mm' or 'm'.

    Returns
    -------
    float
    """
    r_sa = sample_aperture_diameter(input_workspace, unit=unit) / 2.0  # radius
    r_so = source_aperture_diameter(input_workspace, unit=unit) / 2.0  # radius
    l1 = source_sample_distance(input_workspace, unit=unit)
    l2 = sample_detector_distance(input_workspace, unit=unit)

    radius = r_sa + (r_sa + r_so) * (l2 / l1)
    logger.notice("Radius calculated from the input workspace = {:.2} mm".format(radius * 1e3))
    return radius
