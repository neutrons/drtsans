from mantid.simpleapi import (mtd, LoadNexusMonitors)
from drtsans.settings import amend_config
from drtsans.samplelogs import SampleLogs
from drtsans.load import load_events as generic_load_events, sum_data, load_and_split
from drtsans.beam_finder import center_detector, find_beam_center
from drtsans.tof.eqsans.geometry import source_monitor_distance
from drtsans.tof.eqsans.correct_frame import (correct_detector_frame, correct_monitor_frame, transform_to_wavelength)
from drtsans.instruments import extract_run_number, instrument_enum_name
import os

__all__ = ['load_events', 'load_events_monitor', 'sum_data', 'load_events_and_histogram', 'load_and_split']


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
    correct_monitor_frame(output_workspace)
    return mtd[output_workspace]


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

    # Correct TOF of detector
    correct_detector_frame(output_workspace, path_to_pixel=path_to_pixel)

    return mtd[output_workspace]


def load_events_and_histogram(run, detector_offset=0., sample_offset=0., path_to_pixel=True,
                              data_dir=None, output_workspace=None, output_suffix='',
                              center_x=None, center_y=None, mask=None,
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
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace

    """

    # If needed convert comma separated string list of workspaces in list of strings
    if isinstance(run, str):
        run = [r.strip() for r in run.split(',')]

    # if only one run just load and transform to wavelength adn return workspace
    if len(run) == 1:
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

        ws = transform_to_wavelength(ws)
        return ws
    else:
        instrument_unique_name = instrument_enum_name(run[0])  # determine which SANS instrument

        # create default name for output workspace, uses all input
        if (output_workspace is None) or (not output_workspace) or (output_workspace == 'None'):
            output_workspace = '{}_{}{}'.format(instrument_unique_name,
                                                '_'.join(str(extract_run_number(r)) for r in run),
                                                output_suffix)

        # list of worksapce to sum
        temp_workspaces = []

        # load and transform each workspace in turn
        for n, r in enumerate(run):
            temp_workspace_name = '__tmp_ws_{}'.format(n)
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
            transform_to_wavelength(temp_workspace_name)
            temp_workspaces.append(temp_workspace_name)

        # Sum temporary loaded workspaces
        ws = sum_data(temp_workspaces,
                      output_workspace=output_workspace)

        # Remove temporary wokspace
        for ws_name in temp_workspaces:
            if mtd.doesExist(ws_name):
                mtd.remove(ws_name)

        return ws
