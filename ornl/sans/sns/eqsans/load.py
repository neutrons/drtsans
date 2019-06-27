from mantid.simpleapi import LoadEventNexus
from ornl.settings import (optional_output_workspace, amend_config,
                           unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (source_sample_distance,
                                sample_detector_distance)
from ornl.sans.sns.eqsans.geometry import translate_detector_z
from ornl.sans.sns.eqsans.correct_frame import correct_detector_frame

__all__ = ['load_events']


@optional_output_workspace
def load_events(run, output_workspace=None, **levn):
    r"""
    Load events

    Note: Detector is translated along the Z-axis by the value
    specified in keyword "detectorZ" of the logs

    This function is a simplified call to mantid algorithm LoadEventNexus

    Parameters
    ----------
    run: int, str
        Examples: 55555 or EQSANS_55555 or file path
    output_workspace: str
        Name of the output workspace. If none, then use the name of the
        variable
    levn: dict
        Additional positional arguments for LoadEventNexus

    Returns
    -------
    EventWorkspace
        Reference to the events workspace
    """
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        _ws = LoadEventNexus(Filename=str(run),
                             OutputWorkspace=uwd(), **levn)
    #
    # Correct geometry
    #
    # issue #72 LoadInstrument here, loading future EQSANS IDF
    # that moves the detector according to the logs
    translate_detector_z(_ws)  # search logs and translate
    sl = SampleLogs(_ws)
    sl.insert('source-sample-distance', source_sample_distance(_ws), unit='mm')
    sl.insert('sample-detector-distance',
              sample_detector_distance(_ws), unit='mm')
    #
    # Correct TOF
    #
    correct_detector_frame(_ws)
    return _ws
