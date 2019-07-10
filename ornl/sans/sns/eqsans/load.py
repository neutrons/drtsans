
from mantid.simpleapi import (mtd, LoadEventNexus, CloneWorkspace)
from ornl.settings import (amend_config, unique_workspace_name)
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (source_sample_distance,
                                sample_detector_distance)
from ornl.sans.sns.eqsans.geometry import (translate_detector_z,
                                           translate_detector_by_z,
                                           translate_sample_by_z)
from ornl.sans.sns.eqsans.correct_frame import correct_detector_frame
import os

__all__ = ['load_events']


def load_events(run, detector_offset=0., sample_offset=0.,
                output_workspace=None, **kwargs):
    r"""
    Load events with initial corrections for geometry and time-of-flight

    Note: Detector is translated along the Z-axis by the value
    specified in keyword "detectorZ" of the logs

    This function contains a call to mantid algorithm LoadEventNexus.

    Parameters
    ----------
    run: int, str
        Examples: 55555 or EQSANS_55555 or file path.
    detector_offset: float
        Additional translation of the detector along the Z-axis, in mm.
    sample_offset: float
        Additional translation of the sample, in mm. The sample flange remains
        at the origin of coordinates.
    kwargs: dict
        Additional positional arguments for LoadEventNexus.

    Returns
    -------
    EventWorkspace
        Reference to the events workspace
    """
    if output_workspace is None:
        if isinstance(run, str):
            output_workspace = os.path.split(run)[-1]
            output_workspace = '_'.join(output_workspace.split('_')[:2])
            output_workspace = output_workspace.split('.')[0]
        else:
            output_workspace = unique_workspace_name(suffix=str(run))

    if isinstance(run, int) or isinstance(run, str):
        with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
            LoadEventNexus(Filename=str(run),
                           OutputWorkspace=output_workspace, **kwargs)
    else:
        CloneWorkspace(run, OutputWorkspace=output_workspace)
    #
    # Correct geometry
    #
    # issue #72 LoadInstrument here, loading future EQSANS IDF
    # that moves the detector according to the logs
    translate_detector_z(output_workspace)  # search logs and translate
    translate_detector_by_z(output_workspace, 1e-3 * detector_offset)
    translate_sample_by_z(output_workspace, 1e-3 * sample_offset)

    sl = SampleLogs(output_workspace)
    sl.insert('source-sample-distance',
              source_sample_distance(output_workspace, search_logs=False),
              unit='mm')
    sl.insert('sample-detector-distance',
              sample_detector_distance(output_workspace, search_logs=False),
              unit='mm')
    #
    # Correct TOF
    #
    correct_detector_frame(output_workspace)
    return mtd[output_workspace]
