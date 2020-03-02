from mantid.api import MatrixWorkspace
from mantid.geometry import Instrument
from mantid.kernel import logger
# https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.html
from mantid.simpleapi import mtd, MoveInstrumentComponent
import numpy as np

from drtsans.samplelogs import SampleLogs
from drtsans.instruments import InstrumentEnumName, instrument_enum_name
from collections import defaultdict

__all__ = ['beam_radius', 'sample_aperture_diameter', 'source_aperture_diameter', 'translate_sample_by_z',
           'translate_detector_by_z']
detector_z_log = 'detectorZ'


def panel_names(input_query):
    r"""
    List of names for the double-panel detector arrays (e.g., 'detector1', 'wing_detector')

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    list
    """
    detector_names = ['detector1']
    if instrument_enum_name(input_query) == InstrumentEnumName.BIOSANS:
        detector_names.append('wing_detector')
    return detector_names


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
        MatrixWorkspace, workspace name

    Returns
    -------
    Mantid::Instrument
        Instrument object
    """
    def from_ws(ws):
        return ws.getInstrument()

    def from_string(s):
        if s in mtd:
            return get_instrument(mtd[s])

    dispatch = {MatrixWorkspace: from_ws, str: from_string}
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
                    'source_aperture_sample_distance',
                    'source_aperture_sample_aperture_distance')
        if log_key is not None:
            log_keys = (log_key)
        sample_logs = SampleLogs(source)
        try:
            lk = set(log_keys).intersection(set(sample_logs.keys())).pop()
            lk_value = float(sample_logs.single_value(lk))
            # Default unit of lk is mm unless "m" specified
            return lk_value * m2units[unit] if sample_logs[lk].units == 'm' else lk_value * mm2units[unit]
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
        sample_logs = SampleLogs(source)
        try:
            lk = set(log_keys).intersection(set(sample_logs.keys())).pop()
            lk_value = float(sample_logs.single_value(lk))
            # Default unit of lk is mm unless "m" specified
            return lk_value * m2units[unit] if sample_logs[lk].units == 'm' else lk_value * mm2units[unit]
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

    # There are runs for GPSANS and BIOSANS containing log entry "sample_aperture_radius" storing the diameter!
    if 'sample_aperture_radius' in SampleLogs(input_workspace).keys():
        run_limit = {InstrumentEnumName.GPSANS: 7460,
                     InstrumentEnumName.BIOSANS: 1791}.get(instrument_enum_name(input_workspace), 0)
        if int(SampleLogs(input_workspace).run_number.value) < run_limit:
            diameter = SampleLogs(input_workspace).single_value('sample_aperture_radius')

    if diameter is None:
        raise RuntimeError('Unable to retrieve the sample aperture diameter from the logs')

    # If the diameter was found using the additional logs, then insert a log for the diameter under key
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
    sample_logs = SampleLogs(input_workspace)
    diameter = None

    try:
        diameter = sample_logs.single_value('source_aperture_diameter')
    except RuntimeError:
        # There are runs for GPSANS and BIOSANS containing log entry "source_aperture_radius" storing the diameter!
        if 'source_aperture_radius' in sample_logs.keys():
            run_limit = {InstrumentEnumName.GPSANS: 7460,
                         InstrumentEnumName.BIOSANS: 1791}.get(instrument_enum_name(input_workspace), 0)
            if int(SampleLogs(input_workspace).run_number.value) < run_limit:
                diameter = sample_logs.single_value('source_aperture_radius')

    if diameter is None:
        raise ValueError('Unable to retrieve the source aperture diameter from the logs')

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


def translate_source_by_z(input_workspace, z=None, relative=False):
    r"""
      Adjust the Z-coordinate of the source.


      Parameters
      ----------
      input_workspace: ~mantid.api.MatrixWorkspace
          Input workspace containing instrument file
      z: float
          Translation to be applied, in units of meters. If :py:obj:`None`, the quantity stored in the logs
           is used, unless the source has already been translated by this
          quantity.
      relative: bool
          If :py:obj:`True`, add to the current z-coordinate. If :py:obj:`False`, substitute
          the current z-coordinate with the new value.
      """
    if z is None:
        sample_logs = SampleLogs(input_workspace)
        # If detector_z_log exists in the sample logs, use it
        source_z_log = None
        for logname in ['source-sample-distance', 'source_aperture_sample_aperture_distance']:
            if logname in sample_logs:
                source_z_log = logname
                break

        if source_z_log is not None:
            factor = 1.0 if sample_logs[source_z_log].units == 'm' else 1e-3
            distance_from_log = factor * sample_logs.single_value(source_z_log)  # assumed in millimeters
            # Has the detector already been translated by this quantity?
            for source_name in ('moderator', 'source'):
                moderator = get_instrument(input_workspace).getComponentByName(source_name)
                if moderator is not None:
                    _, _, current_z = moderator.getPos()
                    if abs(distance_from_log - abs(current_z)) > 1e-03:  # differ by more than one millimeter
                        z = -distance_from_log
                    break

    if z is not None:
        if (not relative) or (z != 0.):
            for source_name in ('moderator', 'source'):
                if get_instrument(input_workspace).getComponentByName(source_name) is not None:
                    MoveInstrumentComponent(Workspace=input_workspace, Z=z, ComponentName=source_name,
                                            RelativePosition=relative)
                    break


def translate_sample_by_z(workspace, z):
    r"""
    Shift the position of the sample by the desired amount

    Parameters
    ----------
    workspace: ~mantid.api.MatrixWorkspace
        Input workspace containing instrument file
    z: float
        Translation to be applied in meters. Positive values are downstream.
    """
    # only move if the value is non-zero
    if z != 0.:
        MoveInstrumentComponent(Workspace=str(workspace), Z=z,
                                ComponentName='sample-position',
                                RelativePosition=True)

    # update the appropriate log
    sample_logs = SampleLogs(workspace)
    logname_to_set = 'source-sample-distance'  # default
    # look for name of the log/property to update
    for logname in ['source-sample-distance', 'source_aperture_sample_aperture_distance']:
        if logname in sample_logs:
            logname_to_set = logname
            break

    sample_logs.insert(logname_to_set, source_sample_distance(workspace, search_logs=False, unit='mm'),
                       unit='mm')


def translate_detector_by_z(input_workspace, z=None, relative=True):
    r"""
    Adjust the Z-coordinate of the detector.


    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace
        Input workspace containing instrument file
    z: float
        Translation to be applied, in units of meters. If :py:obj:`None`, the quantity stored in log_key
        ~drtsans.geometry.detector_z_log is used, unless the detector has already been translated by this
        quantity.
    relative: bool
        If :py:obj:`True`, add to the current z-coordinate. If :py:obj:`False`, substitute
        the current z-coordinate with the new value.
    """
    update_log = False
    if z is None:
        sample_logs = SampleLogs(input_workspace)
        # If detector_z_log exists in the sample logs, use it
        if detector_z_log in sample_logs:
            translation_from_log = 1e-3 * sample_logs.single_value(detector_z_log)  # assumed in millimeters
            # Has the detector already been translated by this quantity?
            main_detector_array = detector_name(input_workspace)
            _, _, current_z = get_instrument(input_workspace).getComponentByName(main_detector_array).getPos()
            if abs(translation_from_log - current_z) > 1e-03:  # differ by more than one millimeter
                z = translation_from_log

    if z is not None:
        update_log = True
        if (not relative) or (z != 0.):
            MoveInstrumentComponent(Workspace=input_workspace, Z=z, ComponentName=detector_name(input_workspace),
                                    RelativePosition=relative)

    # update the appropriate log
    if update_log:
        sample_logs = SampleLogs(input_workspace)
        sample_logs.insert('sample-detector-distance', sample_detector_distance(input_workspace, search_logs=False),
                           unit='mm')
