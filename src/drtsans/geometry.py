# local imports
from drtsans.samplelogs import SampleLogs
from drtsans.settings import unpack_v3d, namedtuplefy
from drtsans.instruments import InstrumentEnumName, instrument_enum_name, instrument_standard_name
from collections import defaultdict

# third-party imports
from mantid.api import MatrixWorkspace
from mantid.api import Workspace as MantidWorkspace
from mantid.dataobjects import EventWorkspace
from mantid.geometry import Instrument
from mantid.kernel import logger
from mantid.simpleapi import mtd, MoveInstrumentComponent
import numpy as np

# standard imports
from typing import Tuple, Union


__all__ = [
    "sample_aperture_diameter",
    "source_aperture_diameter",
    "translate_sample_by_z",
    "translate_detector_by_z",
    "source_sample_distance",
    "sample_detector_distance",
    "search_sample_detector_distance_meta_name",
    "search_source_sample_distance_meta_name",
]
detector_z_log = "detectorZ"


def panel_names(input_query):
    r"""
    List of names for the double-panel detector arrays (e.g., 'detector1', 'wing_detector')

    Parameters
    ----------
    input_query: str,  ~mantid.api.Workspace
        string representing a valid instrument name, the name of a Mantid workspace,
        or the handle to a Mantid workspace

    Returns
    -------
    list
    """
    # latest panels for each instrument
    panels_for_instrument = {
        "BIOSANS": ["detector1", "wing_detector", "midrange_detector"],
        "GPSANS": ["detector1"],
        "EQSANS": ["detector1"],
    }
    instrument_name = instrument_standard_name(str(input_query))
    panels = panels_for_instrument[instrument_name]
    # the only troublemaker is BIOSANS, which may be missing the midrange detector
    if instrument_name == "BIOSANS" and str(input_query) in mtd:
        if mtd[str(input_query)].getInstrument().getComponentByName("midrange_detector") is None:
            return ["detector1", "wing_detector"]
    return panels


def main_detector_name(ipt):
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
    inst_to_det = defaultdict(lambda: "detector1")
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
    return get_instrument(source).getComponentByName(main_detector_name(source))


def bank_workspace_index_range(input_workspace, component=""):
    """
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
    """
    detector_ids = bank_detector_ids(input_workspace, component, None)
    detector_id_first = detector_ids.min()

    input_workspace = mtd[str(input_workspace)]
    first = None
    for i in range(input_workspace.getNumberHistograms()):
        ids = input_workspace.getSpectrum(i).getDetectorIDs()
        if len(ids) > 1:
            raise RuntimeError("do not know how to work with more than one detector per spectrum ({})".format(ids))
        if ids[0] == detector_id_first:
            first = i
            break
    if first is None:
        raise RuntimeError("something meaningful goes here")
    else:
        return (first, first + detector_ids.size)


def bank_detector_ids(input_workspace, component="", masked=None):
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


def pixel_centers(input_workspace, indexes, shape=None):
    r"""
    Coordinates for the center of one or more pixel detectors. It is assumed that all pixels have the same shape.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    indexes: int, list
        one index or a list of component info indexes, or detector info indexes, but not workspace indexes or
        spectrum info indexes.
    shape: str
        vestigial argument, retained for backward compatiblity

    Returns
    -------
    numpy.ndarray
        Coordinates for each pixel detector
    """
    # NOTE: All pixel positions are defined within IDF with respect to the lab reference frame.
    #       The current IDF for CG2 defines all pixels by their center, which is directly readin
    #       by Mantid upon data loading.
    component_info = mtd[str(input_workspace)].componentInfo()
    indexes = (
        [
            indexes,
        ]
        if isinstance(indexes, int)
        else indexes
    )
    positions = np.array([unpack_v3d(component_info.position, i) for i in indexes])
    return positions


@namedtuplefy
def logged_smearing_pixel_size(input_workspace):
    """Find pixel size (X and Y) from the metadata within a workspace

    Parameters
    ----------
    workspace:  str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        workspace for detector size

    Returns
    -------
    namedtuple
        Fields of the named tuple are, in this order
        - width, float or :py:obj:`None` if no pixel width is found in the logs
        - height, float or :py:obj:`None` if no pixel height is found in the logs
    """
    # Get workspace
    workspace = mtd[str(input_workspace)]

    # Get sample logs
    sample_logs = SampleLogs(workspace)

    if "smearingPixelSizeX" in sample_logs.keys() and "smearingPixelSizeY" in sample_logs.keys():
        smearing_pixel_size_x = sample_logs["smearingPixelSizeX"].value
        smearing_pixel_size_y = sample_logs["smearingPixelSizeY"].value
    else:
        smearing_pixel_size_x, smearing_pixel_size_y = None, None

    return {"width": smearing_pixel_size_x, "height": smearing_pixel_size_y}


@namedtuplefy
def nominal_pixel_size(input_workspace):
    """Find pixel size (X and Y) from the instrument geometry disregarding pixel calibrations like
    barscan and tube-width.

    We use the pixel dimensions of the first detector that is not a monitor detector.

    Parameters
    ----------
    input_workspace:  str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        workspace for detector size

    Returns
    -------
    namedtuple
        Fields of the named tuple are 'width' and 'height', in this order.
    """
    workspace = mtd[str(input_workspace)]  # handle to Mantid Workspace object
    workspace_index = 0
    detectorInfo = workspace.detectorInfo()
    while True:
        if detectorInfo.isMonitor(workspace_index) is False:
            detector = workspace.getDetector(workspace_index)
            detector_shape = detector.shape().getBoundingBox().width()  # (X, Y, Z) values
            return {"width": detector_shape.X(), "height": detector_shape.Y()}
        workspace_index += 1


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


def source_sample_distance(source, unit="mm", log_key=None, search_logs=True):
    r"""
    Report the distance (always positive!) between source and sample aperture.

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
        # log_keys = ('source-sample-distance', 'source_sample-distance',
        #             'source_sample_distance', 'sample-source-distance',
        #             'sample_source-distance', 'sample_source_distance',
        #             'source_aperture_sample_distance',
        #             'source_aperture_sample_aperture_distance')
        # if log_key is not None:
        #     log_keys = (log_key)
        # sample_logs = SampleLogs(source)
        # try:
        #     lk = set(log_keys).intersection(set(sample_logs.keys())).pop()
        #     lk_value = float(sample_logs.single_value(lk))
        #     # Default unit of lk is mm unless "m" specified
        #     return lk_value * m2units[unit] if sample_logs[lk].units == 'm' else lk_value * mm2units[unit]
        # except KeyError:
        #     pass
        # Search the logs for the distance
        meta_info_list = search_source_sample_distance_meta_name(source, log_key)
        if len(meta_info_list) == 0:
            # No meta data found: use instrument geometry to calculate distance
            pass
        else:
            # Retrieve source-sample distance from meta data and convert to correct unit
            # In case that there is more than 1 log that is found, it is assumed that one of them is the native
            # motor position and rest of them are aliases
            # Take the first one shall work
            meta_data_name, meta_distance, meta_data_unit = meta_info_list[0]
            distance = meta_distance * m2units[unit] if meta_data_unit == "m" else meta_distance * mm2units[unit]
            return distance

    # Calculate the distance using the instrument definition file
    instrument = get_instrument(source)
    sample = instrument.getSample()

    return abs(sample.getDistance(instrument.getSource())) * m2units[unit]


def search_source_sample_distance_meta_name(source, specified_meta_name):
    """Search meta data (sample logs) for source-sample distance

    Parameters
    ----------
    source : PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number.
    specified_meta_name : str, None
        Only search for the given string in the source's meta data (logs). Do not use default log keys

    Returns
    -------
    ~list
        item = (str, float, str)
        meta data name, sample detector distance value, unit

    """
    # Allowed meta data name for source-sample distance including aliases
    # Currently EPICS uses 'source_aperture_sample_aperture_distance'
    log_key_order = [
        "source_aperture_sample_aperture_distance",
        "source_aperture_sample_distance",
        "source-sample-distance",
        "source_sample-distance",
        "source_sample_distance",
        "sample-source-distance",
        "sample_source-distance",
        "sample_source_distance",
    ]
    log_keys = set(log_key_order)
    meta_data_list = _search_meta_data(source, log_keys, specified_meta_name)

    # return the index of meta_data_name in the log_key_order list
    def get_log_index(meta_data):
        return log_key_order.index(meta_data[0])

    # the meta_data are sorted based on log_key_order items' order
    sorted_meta_data_list = sorted(meta_data_list, key=get_log_index)
    return sorted_meta_data_list


def sample_detector_distance(
    source,
    unit: str = "mm",
    log_key: Union[str, None] = None,
    search_logs: bool = True,
    forbid_calculation: bool = False,
) -> float:
    r"""
    Return the distance from the sample to the main detector bank plane

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
    forbid_calculation: bool
        Flag to raise an exception if it is required to get SDD from meta data but no associated meta data is found

    Returns
    -------
    float
        distance between sample and detector, in selected units
    """
    m2units = dict(mm=1e3, m=1.0)
    mm2units = dict(mm=1.0, m=1e-3)

    if search_logs is True:
        # Search the logs for the distance
        meta_info_list = search_sample_detector_distance_meta_name(source, log_key)
        if len(meta_info_list) == 0:
            # No meta data found: Use instrument information to get distance
            if forbid_calculation:
                raise RuntimeError("Unable to find any meta data related to SDD")
        else:
            # Calculate from log value considering unit
            # In case there are more than 1 log is found, it is assumed that all of them shall have the same
            # value, i.e., some of them are alias
            meta_data_name, meta_data_value, meta_data_unit = meta_info_list[0]
            distance = meta_data_value * m2units[unit] if meta_data_unit == "m" else meta_data_value * mm2units[unit]
            return distance

    # Calculate the distance using the instrument definition file
    instrument = get_instrument(source)
    det = instrument.getComponentByName(main_detector_name(source))

    # Mantid instrument uses 'meter' only
    det_pos = det.getPos()
    sample_pos = instrument.getSample().getPos()
    sample_det_plane_distance = abs(det_pos.Z() - sample_pos.Z())

    return sample_det_plane_distance * m2units[unit]


def search_sample_detector_distance_meta_name(source, specified_meta_name):
    """Search meta data (sample logs) for sample detector distance

    Parameters
    ----------
    source : PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number.
    specified_meta_name : str, None
        Only search for the given string in the source's meta data (logs). Do not use default log keys

    Returns
    -------
    ~list
        item = (str, float, str)
        meta data name, sample detector distance value, unit

    """
    log_keys = {
        "detector-sample-distance",
        "detector_sample-distance",
        "detector_sample_distance",
        "sample-detector-distance",
        "sample_detector-distance",
        "sample_detector_distance",
    }  # latest one

    return _search_meta_data(source, log_keys, specified_meta_name)


def _search_meta_data(source, default_search_set, specified_meta_name):
    """

    Parameters
    ----------
    source : PyObject
        Instrument object, MatrixWorkspace, workspace name, file name,
        run number.
    default_search_set: ~set
        Set of strings of possible meta data name in workspace's run properties
    specified_meta_name: str
        User-specified meta data name overriding default_search_set
    Returns
    -------
    ~list
        item = (str, float, str)
        meta data name, sample detector distance value, unit

    """
    # Determine the possible meta data names to search for
    if specified_meta_name is None:
        # Possible meta data name
        log_keys = default_search_set
    else:
        # User specified
        log_keys = {specified_meta_name}
    # Get sample logs
    sample_logs = SampleLogs(source)

    # Intersection:
    found_log_names = sorted(list(log_keys.intersection((sample_logs.keys()))))

    # Decide log name
    meta_list = list()
    for meta_name in found_log_names:
        lk_value = float(sample_logs.single_value(meta_name))
        distance_unit = sample_logs[meta_name].units
        meta_list.append((meta_name, lk_value, distance_unit))

    return meta_list


def source_detector_distance(source, unit="mm", search_logs=True):
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


def sample_aperture_diameter(input_workspace, unit="mm"):
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
    log_keys = ["sample_aperture_diameter"]

    try:
        additional_log_keys = {
            InstrumentEnumName.EQSANS: ["beamslit4"],
            InstrumentEnumName.GPSANS: [],
            InstrumentEnumName.BIOSANS: [],
        }
        log_keys += additional_log_keys[instrument_enum_name(input_workspace)]
    except KeyError:
        # In case the instrument name (test instrument) not in EQ, GP and BIO-SANS
        pass

    sample_logs = SampleLogs(input_workspace)
    diameter = None

    for log_key in log_keys:
        if log_key in sample_logs.keys():
            diameter = sample_logs.single_value(log_key)
            break

    # There are runs for GPSANS and BIOSANS containing log entry "sample_aperture_radius" storing the diameter!
    if "sample_aperture_radius" in SampleLogs(input_workspace).keys():
        run_limit = {
            InstrumentEnumName.GPSANS: 7460,
            InstrumentEnumName.BIOSANS: 1791,
        }.get(instrument_enum_name(input_workspace), 0)
        if int(SampleLogs(input_workspace).run_number.value) < run_limit:
            diameter = SampleLogs(input_workspace).single_value("sample_aperture_radius")

    if diameter is None:
        raise RuntimeError("Unable to retrieve the sample aperture diameter from the logs")

    # If the diameter was found using the additional logs, then insert a log for the diameter under key
    # "sample_aperture_diameter"
    if "sample_aperture_diameter" not in sample_logs.keys():
        sample_logs.insert("sample_aperture_diameter", diameter, unit="mm")

    return diameter if unit == "mm" else diameter / 1.0e3


def source_aperture_diameter(input_workspace, unit="mm"):
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
        diameter = sample_logs.single_value("source_aperture_diameter")
    except RuntimeError:
        # There are runs for GPSANS and BIOSANS containing log entry "source_aperture_radius" storing the diameter!
        if "source_aperture_radius" in sample_logs.keys():
            run_limit = {
                InstrumentEnumName.GPSANS: 7460,
                InstrumentEnumName.BIOSANS: 1791,
            }.get(instrument_enum_name(input_workspace), 0)
            if int(SampleLogs(input_workspace).run_number.value) < run_limit:
                diameter = sample_logs.single_value("source_aperture_radius")

    if diameter is None:
        raise ValueError("Unable to retrieve the source aperture diameter from the logs")

    return diameter if unit == "mm" else diameter / 1.0e3


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
        for logname in [
            "source-sample-distance",
            "source_aperture_sample_aperture_distance",
        ]:
            if logname in sample_logs:
                source_z_log = logname
                break

        if source_z_log is not None:
            factor = 1.0 if sample_logs[source_z_log].units == "m" else 1e-3
            distance_from_log = factor * sample_logs.single_value(source_z_log)  # assumed in millimeters
            # Has the detector already been translated by this quantity?
            for source_name in ("moderator", "source"):
                moderator = get_instrument(input_workspace).getComponentByName(source_name)
                if moderator is not None:
                    _, _, current_z = moderator.getPos()
                    if abs(distance_from_log - abs(current_z)) > 1e-03:  # differ by more than one millimeter
                        z = -distance_from_log
                    break

    if z is not None:
        if (not relative) or (z != 0.0):
            for source_name in ("moderator", "source"):
                if get_instrument(input_workspace).getComponentByName(source_name) is not None:
                    MoveInstrumentComponent(
                        Workspace=input_workspace,
                        Z=z,
                        ComponentName=source_name,
                        RelativePosition=relative,
                    )
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
    if z != 0.0:
        ws_name = str(workspace)
        MoveInstrumentComponent(
            Workspace=str(workspace),
            Z=z,
            ComponentName="sample-position",
            RelativePosition=True,
        )
        workspace = mtd[ws_name]
        logger.debug(
            "Instrument sample position is moved to {}".format(workspace.getInstrument().getSample().getPos())
        )

    # update the appropriate log
    # 'source_aperture_sample_aperture_distance' is not coupled with sample/source distance. Thus
    # it won't be updated
    # FIXME - This is a technical debt because NO instrument or other drt-sans code actually does use this value
    sample_logs = SampleLogs(workspace)
    logname_to_set = "source-sample-distance"  # default
    sample_logs.insert(
        logname_to_set,
        source_sample_distance(workspace, search_logs=False, unit="mm"),
        unit="mm",
    )
    # FIXME - sample-detector distance shall be updated too


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
            main_detector_array = main_detector_name(input_workspace)
            _, _, current_z = get_instrument(input_workspace).getComponentByName(main_detector_array).getPos()
            if abs(translation_from_log - current_z) > 1e-03:  # differ by more than one millimeter
                z = translation_from_log

    if z is not None:
        update_log = True
        if (not relative) or (z != 0.0):
            logger.debug(
                "Moving detector along Z = {}  is relative = {} to component {}".format(
                    z, relative, main_detector_name(input_workspace)
                )
            )

            MoveInstrumentComponent(
                Workspace=input_workspace,
                Z=z,
                ComponentName=main_detector_name(input_workspace),
                RelativePosition=relative,
            )

    # update the appropriate log
    if update_log:
        sample_logs = SampleLogs(input_workspace)
        sample_logs.insert(
            "sample-detector-distance",
            sample_detector_distance(input_workspace, search_logs=False),
            unit="mm",
        )


#################################################################################
#   Functionality to help simulate events or intensities for a component panel
#################################################################################

PIXELS_IN_TUBE: int = 256


def get_position_south_detector(input_workspace: Union[str, MantidWorkspace, EventWorkspace]) -> float:
    r"""
    Get the downstream position of the south detector.

    Parameters
    ----------
    Name or handle to a Mantid Workspace.

    Returns
    -------
    Position of the south detector from the origin, in meters
    """
    workspace = mtd[str(input_workspace)]
    mantid_instrument = workspace.getInstrument()
    south_detector = mantid_instrument.getComponentByName("detector1")
    return south_detector.getPos().Z()


def get_curvature_radius(input_workspace: Union[str, MantidWorkspace], component: str) -> float:
    r"""
    Get the curvature radius of the given component.

    Currently, only the BIOSANS instrument is supported since it's the only one with curved detectors.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the curved double-panels in the instrument (e.g. "wing_detector", "midrange_detector")

    Returns
    -------
    Nominal radius of curvature of the given component, in meters
    """
    instrument = instrument_standard_name(input_workspace)
    radii = {
        "BIOSANS": {
            "wing_detector": 1.1633,
            "midrange_detector": 4.0000,
        }
    }
    return radii[instrument][component]


def spectrum_info_ranges(
    input_workspace: Union[str, MantidWorkspace, EventWorkspace], component: str
) -> Tuple[int, int]:
    r"""
    Get the spectra range (as workspace indexes) for the given component. This assumes monitors have no
    spectra associated to them.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")

    Returns
    -------
    First and next-to-last workspace indexes for the given component
    """
    ranges = {
        "BIOSANS": {
            "detector1": (0, 49152),  # 49152 is the first pixel in the next component
            "wing_detector": (49152, 90112),  # 90112 is the first pixel in the next component
            "midrange_detector": (90112, 106496),
        },
        "EQSANS": {
            "detector1": (0, 49152),  # 49152 is 1 + last_pixel_index
        },
        "GPSANS": {
            "detector1": (0, 49152),  # 49152 is 1 + last_pixel_index
        },
    }
    instrument = instrument_standard_name(input_workspace)
    return ranges[instrument][component]


def get_pixel_masks(input_workspace: Union[str, MantidWorkspace], component: str) -> np.ndarray:
    r"""
    Get the pixel masks for the given component.

    Parameters
    ----------
    input_workspace
            Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")

    Returns
    -------
    Array of True/False values, where True indicates that the pixel is masked.
    """
    spectrum_info = mtd[str(input_workspace)].spectrumInfo()
    first, next_to_last = spectrum_info_ranges(input_workspace, component)
    return np.array([spectrum_info.isMasked(idx) for idx in range(first, next_to_last)])


def get_pixel_distances(input_workspace: Union[str, MantidWorkspace], component: str) -> np.ndarray:
    r"""
    Get the distances from the sample to each of the pixel detectors in the given component, in meters.

    Parameters
    ----------
    input_workspace
        The workspace to calculate the solid angles for.

    component
        The component to calculate the solid angles for. (e.g "midrange_detector", "detector1").

    Returns
    -------
    np.ndarray
    """
    spectrum_info = mtd[str(input_workspace)].spectrumInfo()
    first, next_to_last = spectrum_info_ranges(input_workspace, component)
    return np.array([spectrum_info.l2(idx) for idx in range(first, next_to_last)])


def get_xyz(input_workspace: Union[str, MantidWorkspace], component: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Get the X, Y, and Z coordinates of the pixels in the given component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")
    Returns
    -------
    Tuple of all X, Y, and Z coordinates of the pixels in the given component
    """
    spectrum_info = mtd[str(input_workspace)].spectrumInfo()
    first, next_to_last = spectrum_info_ranges(input_workspace, component)
    x = np.array([spectrum_info.position(idx).X() for idx in range(first, next_to_last)])
    y = np.array([spectrum_info.position(idx).Y() for idx in range(first, next_to_last)])
    z = np.array([spectrum_info.position(idx).Z() for idx in range(first, next_to_last)])

    return x, y, z


def get_xy(input_workspace: Union[str, MantidWorkspace], component: str) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Get the X and Y coordinates of the pixels in the given component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")
    Returns
    -------
    Tuple of all X and Y coordinates of the pixels in the given component
    """
    x, y, z = get_xyz(input_workspace, component)
    return x, y


_instrument_solid_angles_cache = {
    "BIOSANS": {"midrange_detector": None, "wing_detector": None, "detector1": None},
    "EQSANS": {"detector1": None},
    "GPSANS": {"detector1": None},
}


def get_twothetas(input_workspace: Union[str, MantidWorkspace], component: str, units="degrees") -> np.ndarray:
    r"""
    Get the two-theta angles for each of the pixel detectors in the given component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")
    units
        The units to return the two-theta angles in. Must be one of "degrees" or "radians".

    Returns
    -------
    Array of two-theta angles for each of the pixel detectors in the given component.
    """
    if units.lower() not in ["degrees", "radians"]:
        raise ValueError("units must be 'degrees' or 'radians'")
    spectrum_info = mtd[str(input_workspace)].spectrumInfo()
    first, next_to_last = spectrum_info_ranges(input_workspace, component)
    two_thetas = np.array([spectrum_info.twoTheta(idx) for idx in range(first, next_to_last)])
    if units.lower() == "degrees":
        return np.degrees(two_thetas)
    else:
        return two_thetas


def get_solid_angles(
    input_workspace: Union[str, MantidWorkspace], component: str, back_panel_attenuation: float = 0.5
) -> np.ndarray:
    r"""
    Calculate (and cache) the solid angles subtended by the detector pixels of the given component.

    Pixel detectors in the Wing and Midrange componets are invariant under rotation of the whole component. Therefore,
    the solid angles are computed once. The solid angles for the South Detector component depend on the position of the
    component, hence solid angles are cached at each position, with a generous resolution of 0.5 m in position.

    Parameters
    ----------
    input_workspace
        The workspace to calculate the solid angles for.
    component
        The component to calculate the solid angles for. Must be one of "midrange_detector", "wing_detector", or
        "detector1".
    back_panel_attenuation
        An approximate constant attenuation factor for pixels in the back panel of the component, since they're
        partially occluded by the pixels in the front panel.

    Returns
    -------
        The solid angles subtended by the detector pixels of the given component.
    """

    def _get_solid_angles():
        r"""The solid angle of a pixel with coordinates (x, y, z) is given by cos(chi) / r^2 = (x^2 + Z^2) / r^3.
        Angle chi is the angle between the position vector of the pixel detector and the horizontal (XZ) plane.

        """
        if component == "detector1":
            # approximation for the South panel: assume it's a curved detector.
            radius = get_position_south_detector(input_workspace)
        else:
            radius = get_curvature_radius(input_workspace, component)

        distances = get_pixel_distances(input_workspace, component)

        # pixels in the front and back panel of the component alternate every four Helium tubes
        # we first calculate the attenuation in and eightpack (8 Helium tubes)
        pixels_in_a_fourpack = 4 * PIXELS_IN_TUBE
        attenuation = np.concatenate(
            (np.repeat(1.0, pixels_in_a_fourpack), np.repeat(back_panel_attenuation, pixels_in_a_fourpack))
        )
        # repeat the eightpack attenuation for all eightpacks in the component
        pixel_count = len(distances)
        eightpack_count = int(pixel_count / (2 * pixels_in_a_fourpack))
        attenuation = np.tile(attenuation, eightpack_count)
        return attenuation * radius / np.power(distances, 3)

    solid_angles_cache = _instrument_solid_angles_cache[instrument_standard_name(input_workspace)]
    DISTANCE_BIN_WIDTH = 0.5  # meters
    if solid_angles_cache[component] is None:
        if component == "detector1":
            position_bin = int(get_position_south_detector(input_workspace) / DISTANCE_BIN_WIDTH)
            solid_angles_cache[component] = {position_bin: _get_solid_angles()}
            return solid_angles_cache[component][position_bin]
        else:  # wing or midrange
            solid_angles_cache[component] = _get_solid_angles()
            return solid_angles_cache[component]
    else:
        if component == "detector1":
            position_bin = int(get_position_south_detector(input_workspace) / DISTANCE_BIN_WIDTH)
            if position_bin not in solid_angles_cache[component]:
                solid_angles_cache[component][position_bin] = _get_solid_angles()
            return solid_angles_cache[component][position_bin]
        else:  # wing or midrange
            return solid_angles_cache[component]
