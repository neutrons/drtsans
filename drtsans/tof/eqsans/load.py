from mantid.simpleapi import (mtd, LoadNexusMonitors)
from drtsans.settings import amend_config, namedtuplefy
from drtsans.samplelogs import SampleLogs
from drtsans.load import load_events as generic_load_events, sum_data, load_and_split as generic_load_and_split
from drtsans.beam_finder import center_detector, find_beam_center
from drtsans.tof.eqsans.geometry import source_monitor_distance
from drtsans.tof.eqsans.correct_frame import (correct_detector_frame,
                                              correct_monitor_frame, transform_to_wavelength, smash_monitor_spikes,
                                              set_init_uncertainties, correct_tof_offset, correct_emission_time)
from drtsans.instruments import extract_run_number, instrument_enum_name
import os

__all__ = ['load_events', 'load_events_monitor', 'sum_data', 'load_events_and_histogram',
           'load_and_split', 'prepare_monitors']


def load_events_monitor(run, data_dir=None, output_workspace=None):
    r"""
    Load monitor events with initial corrections for geometry
    and time-of-flight

    Parameters
    ----------
    run: int, str
        Examples: ``55555`` or ``EQSANS_55555`` or file path.
    data_dir: str, list
        Additional data search directories
    output_workspace: str
        If not specified it will be ``EQSANS_55555_monitors`` determined from
        the supplied value of ``run``

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    suffix = '_monitors'
    if output_workspace is None:
        if isinstance(run, str):
            output_workspace = os.path.split(run)[-1]
            output_workspace = '_'.join(output_workspace.split('_')[:2])
            output_workspace = output_workspace.split('.')[0] + suffix
        else:
            output_workspace = 'EQSANS_{}{}'.format(run, suffix)

    with amend_config({'default.instrument': 'EQSANS'}, data_dir=data_dir):
        LoadNexusMonitors(Filename=str(run), LoadOnly='Events', OutputWorkspace=output_workspace)

    smd = source_monitor_distance(output_workspace, unit='mm')
    SampleLogs(output_workspace).insert('source-monitor-distance', smd,
                                        unit='mm')

    # Correct TOF offset
    correct_tof_offset(output_workspace)

    # DAS automatically corrects the frame for the monitor only on or after June 1 2019,
    day_stamp = int(SampleLogs(output_workspace).start_time.value[0:10].replace('-', ''))
    if day_stamp < 20190601:
        correct_monitor_frame(output_workspace)

    # Correct TOF for emission time
    correct_emission_time(output_workspace)

    return mtd[output_workspace]


def prepare_monitors(data, bin_width=0.1, output_workspace=None):
    r"""
    Loads monitor counts, correct TOF, and transforms to wavelength.

    Parameters
    ----------
    data: int, str, ~mantid.api.IEventWorkspace
        Run number as int or str, file path, :py:obj:`~mantid.api.IEventWorkspace`
    bin_width: float
        Bin width for the output workspace, in Angstroms.
    output_workspace: str
        Name of the output workspace. If None, then it will be
        ``EQSANS_XXXXX_monitors`` with number XXXXX determined from ``data``.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    w = load_events_monitor(data, output_workspace=output_workspace)
    w = smash_monitor_spikes(w)
    w = transform_to_wavelength(w, bin_width=bin_width)
    w = set_init_uncertainties(w)
    return w


def load_events(run, detector_offset=0., sample_offset=0., path_to_pixel=True,
                data_dir=None, output_workspace=None, output_suffix='', **kwargs):
    r"""
    Load events with initial corrections for geometry and time-of-flight

    Note: Detector is translated along the Z-axis by the value specified in keyword ``detectorZ`` of the logs. The
    final sample to detector distance is detectorZ + detector_offset - sample_offset.

    **Mantid algorithms used:**
    :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`,

    Parameters
    ----------
    run: int, str
        Examples: ``55555`` or ``EQSANS_55555`` or file path.
    detector_offset: float
        Additional translation of the detector along the Z-axis, in mm. Positive
        moves the detector downstream.
    sample_offset: float
        Additional translation of the sample, in mm. The sample flange remains
        at the origin of coordinates. Positive moves the sample downstream.
    path_to_pixel: bool
        When correcting the recorded time of flight of each neutron, use the
        path from the moderator to the detector pixel (`True`) or to the center
        of the detector panel (`False`). The latter for comparison to the
        old data reduction.
    data_dir: str, list
        Additional data search directories
    output_workspace: str
        If not specified it will be ``EQSANS_55555`` determined from the supplied
        value of ``run``
    output_suffix: str
        If the ``output_workspace`` is not specified, this is appended to the automatically generated
        output workspace name.
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    # use the generic functionality to do most of the work
    output_workspace = generic_load_events(run=run, data_dir=data_dir, output_workspace=output_workspace,
                                           output_suffix=output_suffix, detector_offset=detector_offset,
                                           sample_offset=sample_offset)

    # EQSANS specific part benefits from converting workspace to a string
    output_workspace = str(output_workspace)

    # Correct TOF offset
    correct_tof_offset(output_workspace)
    # Correct TOF of detector
    correct_detector_frame(output_workspace, path_to_pixel=path_to_pixel)
    # Correct TOF for emission time
    correct_emission_time(output_workspace)

    return mtd[output_workspace]


@namedtuplefy
def load_events_and_histogram(run, detector_offset=0., sample_offset=0., path_to_pixel=True,
                              data_dir=None, output_workspace=None, output_suffix='',
                              bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                              center_x=None, center_y=None, mask=None, monitors=False,
                              keep_events=True,
                              **kwargs):
    r"""Load events from one or more NeXus files with initial corrections
    for geometry, time-of-flight and beam center. Convert to
    wavelength and sum.

    Note: Detector is translated along the Z-axis by the value
    specified in keyword ``detectorZ`` of the logs. The final sample
    to detector distance is detectorZ + detector_offset -
    sample_offset.

    Parameters
    ----------
    run: list of runs to load
        Examples: ``55555`` or ``EQSANS_55555`` or file path or ``EQSANS_55555, EQSANS_55556``
    detector_offset: float
        Additional translation of the detector along the Z-axis, in mm. Positive
        moves the detector downstream.
    sample_offset: float
        Additional translation of the sample, in mm. The sample flange remains
        at the origin of coordinates. Positive moves the sample downstream.
    path_to_pixel: bool
        When correcting the recorded time of flight of each neutron, use the
        path from the moderator to the detector pixel (`True`) or to the center
        of the detector panel (`False`). The latter for comparison to the
        old data reduction.
    data_dir: str, list
        Additional data search directories
    output_workspace: str
        If not specified it will be ``EQSANS_55555`` determined from the supplied
        value of ``run``
    output_suffix: str
        If the ``output_workspace`` is not specified, this is appended to the automatically generated
        output workspace name.
    bin_width: float
        Bin width for the output workspace, in Angstroms.
    low_tof_clip: float
        Ignore events with a time-of-flight (TOF) smaller than the minimal
        TOF plus this quantity.
    high_tof_clip: float
        Ignore events with a time-of-flight (TOF) bigger than the maximal
        TOF minus this quantity.
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    mask: mask file path, MaskWorkspace, list
        Additional mask to be applied. If `list`, it is a list of
        detector ID's.
    monitors: boolean
        Option to load the monitors as well as the data, if False monitor will be None
    keep_events: bool
        The final histogram will be an EventsWorkspace if True.
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    namedtuple
        Fields of namedtuple
        data: the loaded data
        monitor: the monitor for the data, if monitors==True else None
    """

    # If needed convert comma separated string list of workspaces in list of strings
    if isinstance(run, str):
        run = [r.strip() for r in run.split(',')]

    if output_workspace is not None:
        monitor_workspace = output_workspace + "_monitors"
    else:
        monitor_workspace = None
    # if only one run just load and transform to wavelength and return workspace
    if len(run) == 1:
        if monitors:
            ws_monitors = prepare_monitors(run[0], bin_width, monitor_workspace)
        else:
            ws_monitors = None

        ws = load_events(run=run[0],
                         detector_offset=detector_offset,
                         sample_offset=sample_offset,
                         path_to_pixel=path_to_pixel,
                         data_dir=data_dir,
                         output_workspace=output_workspace,
                         output_suffix=output_suffix,
                         **kwargs)

        if center_x is None or center_y is None:
            center_x, center_y = find_beam_center(ws, mask=mask)
        center_detector(ws, center_x=center_x, center_y=center_y)  # operates in-place

        ws = transform_to_wavelength(ws, bin_width=bin_width,
                                     low_tof_clip=low_tof_clip,
                                     high_tof_clip=high_tof_clip,
                                     keep_events=keep_events)
        ws = set_init_uncertainties(ws)

        return dict(data=ws,
                    monitor=ws_monitors)
    else:
        if keep_events:
            raise NotImplementedError("Cannot merge runs together with keep_events=True.")

        instrument_unique_name = instrument_enum_name(run[0])  # determine which SANS instrument

        # create default name for output workspace, uses all input
        if (output_workspace is None) or (not output_workspace) or (output_workspace == 'None'):
            output_workspace = '{}_{}{}'.format(instrument_unique_name,
                                                '_'.join(str(extract_run_number(r)) for r in run),
                                                output_suffix)

        # list of worksapce to sum
        temp_workspaces = []
        temp_monitors_workspaces = []

        # load and transform each workspace in turn
        for n, r in enumerate(run):
            temp_workspace_name = '__tmp_ws_{}'.format(n)
            temp_monitors_workspace_name = temp_workspace_name+"_monitors"
            if monitors:
                prepare_monitors(r, bin_width, temp_monitors_workspace_name)

            load_events(run=r,
                        detector_offset=detector_offset,
                        sample_offset=sample_offset,
                        path_to_pixel=path_to_pixel,
                        data_dir=data_dir,
                        output_workspace=temp_workspace_name,
                        output_suffix=output_suffix,
                        **kwargs)
            if center_x is None or center_y is None:
                center_x, center_y = find_beam_center(temp_workspace_name, mask=mask)
            center_detector(temp_workspace_name, center_x=center_x, center_y=center_y)  # operates in-place
            transform_to_wavelength(temp_workspace_name,
                                    bin_width=bin_width,
                                    low_tof_clip=low_tof_clip,
                                    high_tof_clip=high_tof_clip,
                                    keep_events=keep_events)
            temp_workspaces.append(temp_workspace_name)

        # Sum temporary loaded monitor workspaces
        if monitors:
            ws_monitors = sum_data(temp_monitors_workspaces,
                                   output_workspace=monitor_workspace)
            # After summing data re-calculate initial uncertainties
            ws_monitors = set_init_uncertainties(ws_monitors)
        else:
            ws_monitors = None

        # Sum temporary loaded workspaces
        ws = sum_data(temp_workspaces,
                      output_workspace=output_workspace)

        # After summing data re-calculate initial uncertainties
        ws = set_init_uncertainties(ws)

        # Remove temporary wokspace
        for ws_name in temp_workspaces:
            if mtd.doesExist(ws_name):
                mtd.remove(ws_name)

        return dict(data=ws,
                    monitor=ws_monitors)


def load_and_split(run, detector_offset=0., sample_offset=0., path_to_pixel=True,
                   data_dir=None, output_workspace=None, overwrite_instrument=True, output_suffix='',
                   bin_width=0.1, low_tof_clip=500, high_tof_clip=2000,
                   center_x=None, center_y=None, mask=None, monitors=False,
                   keep_events=True,
                   time_interval=None, log_name=None, log_value_interval=None,
                   reuse_workspace=False, **kwargs):

    ws = generic_load_and_split(run=run, data_dir=data_dir,
                                output_workspace=output_workspace, overwrite_instrument=overwrite_instrument,
                                output_suffix=output_suffix,
                                detector_offset=detector_offset, sample_offset=sample_offset,
                                time_interval=time_interval, log_name=log_name, log_value_interval=log_value_interval,
                                reuse_workspace=reuse_workspace, monitors=False,
                                instrument_unique_name='EQSANS', **kwargs)

    for _w in ws:
        # Correct TOF offset
        correct_tof_offset(_w)
        # Correct TOF of detector
        correct_detector_frame(_w, path_to_pixel=path_to_pixel)
        # Correct TOF for emission time
        correct_emission_time(_w)
        if center_x is None or center_y is None:
            center_x, center_y = find_beam_center(_w, mask=mask)
        center_detector(_w, center_x=center_x, center_y=center_y)  # operates in-place

        transform_to_wavelength(_w, bin_width=bin_width,
                                low_tof_clip=low_tof_clip,
                                high_tof_clip=high_tof_clip,
                                keep_events=keep_events)
        set_init_uncertainties(_w)

    return ws
