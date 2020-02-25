from mantid.simpleapi import (mtd, LoadNexusMonitors)
from drtsans.settings import amend_config
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import sample_detector_distance
from drtsans.load import load_events as generic_load_events
from drtsans.tof.eqsans.geometry import (translate_detector_z, translate_detector_by_z, source_monitor_distance)
from drtsans.tof.eqsans.correct_frame import (correct_detector_frame, correct_monitor_frame)
import os

__all__ = ['load_events', 'load_events_monitor']


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
        Additional translation of the detector along the Z-axis, in mm.
    sample_offset: float
        Additional translation of the sample, in mm. The sample flange remains
        at the origin of coordinates.
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
                                           output_suffix=output_suffix, sample_offset=-1. * sample_offset)

    # EQSANS specific part benefits from converting workspace to a string
    output_workspace = str(output_workspace)

    # Correct the distances between instrument components
    translate_detector_z(output_workspace)  # search logs and translate if necessary
    translate_detector_by_z(output_workspace, 1e-3 * detector_offset)

    sample_logs = SampleLogs(output_workspace)
    sample_logs.insert('sample-detector-distance', sample_detector_distance(output_workspace, search_logs=False),
                       unit='mm')

    # Correct TOF of detector
    correct_detector_frame(output_workspace, path_to_pixel=path_to_pixel)

    return mtd[output_workspace]
