from mantid.simpleapi import (LoadEventNexus, MoveInstrumentComponent)
from ornl.settings import (optional_output_workspace,
                           unique_workspace_dundername as uwd)
from ornl.sans.sns.eqsans.geometry import translate_detector_z

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
    _ws = LoadEventNexus(Filename=str(run),
                         OutputWorkspace=uwd(), **levn)
    # issue #72 LoadInstrument here, loading future EQSANS IDF
    # that moves the detector according to the logs
    translate_detector_z(_ws)  # search logs and translate
    return _ws
