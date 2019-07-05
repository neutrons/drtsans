
from mantid.simpleapi import (LoadEventNexus, CloneWorkspace)
from ornl.settings import (optional_output_workspace, amend_config,
                           unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (source_sample_distance,
                                sample_detector_distance)
from ornl.sans.sns.eqsans.geometry import (translate_detector_z,
                                           translate_detector_by_z,
                                           translate_sample_by_z)
from ornl.sans.sns.eqsans.correct_frame import correct_detector_frame

__all__ = ['load_events']


@optional_output_workspace
def load_events(run, detector_offset=0., sample_offset=0., **levn):
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
    levn: dict
        Additional positional arguments for LoadEventNexus.

    Returns
    -------
    EventWorkspace
        Reference to the events workspace
    """
    if isinstance(run, int) or isinstance(run, str):
        with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
            _ws = LoadEventNexus(Filename=str(run),
                                 OutputWorkspace=uwd(), **levn)
    else:
        _ws = CloneWorkspace(run, OutputWorkspace=uwd())
    #
    # Correct geometry
    #
    # issue #72 LoadInstrument here, loading future EQSANS IDF
    # that moves the detector according to the logs
    translate_detector_z(_ws)  # search logs and translate
    translate_detector_by_z(_ws, 1e-3 * detector_offset)
    translate_sample_by_z(_ws, 1e-3 * sample_offset)

    sl = SampleLogs(_ws)
    sl.insert('source-sample-distance',
              source_sample_distance(_ws, search_logs=False), unit='mm')
    sl.insert('sample-detector-distance',
              sample_detector_distance(_ws, search_logs=False), unit='mm')
    #
    # Correct TOF
    #
    correct_detector_frame(_ws)
    return _ws
