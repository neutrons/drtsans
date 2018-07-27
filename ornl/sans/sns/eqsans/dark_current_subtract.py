from __future__ import (absolute_import, division, print_function)

from mantid.api import (PythonAlgorithm, AlgorithmFactory, WorkspaceProperty)
from mantid.simpleapi import Scale
from mantid.kernel import Direction, logger


def compute_log_ratio(run1, run2, l):
    '''Compute ratio of data to dark for one log entry
    :param run1: first run object containing logs
    :param run2: second run object containing logs
    :param l: entry log
    :return: ratio
    :except: RuntimeError when log is not present
    '''
    dt = run1.getProperty(l)
    dk = run2.getProperty(l)
    try:
        dt = dt.getStatistics().duration
        dk = dk.getStatistics().duration
    except AttributeError:
        dt = dt.value
        dk = dk.value
    return dt / dk


def duration_ratio(data, dark, logname=None):
    """
    Compute the ratio of data to dark-current durations

    :param data: run object for data
    :param darkcurrent: run object for dark current
    :param logname: entry log containing the duration. If None, duration will
    be tried looking sequentially into log entries duration', 'proton_charge',
    and 'timer.
    :return: duration ratio, or 1.0 if ratio cannot be computed
    """
    def valid_logname(l):
        return data.hasProperty(l) and dark.hasProperty(l)

    not_found = 'Logs could not be found, duration ratio set to 1.0'
    if logname is not None:
        try:
            return compute_log_ratio(data, dark, logname)
        except RuntimeError:
            logger.error(not_found)
    else:
        for l in ('duration', 'proton_charge', 'timer'):
            try:
                return compute_log_ratio(data, dark, l)
            except RuntimeError:
                continue
    return 1.0


def subtract_scaled_dark(data, dark, logname=None):
    """Rescale dark current and then subtract from data.

    Events from the dark current are included in as weighted neutron events.
    The weight for each dark event is -duration_ratio, with duration_ratio
    being the ratio of the data to dark current duration runs.

    The scaling factor should account for the TOF cuts on each side of a frame
    The EQSANSLoad algorithm cuts the beginning and end of the TOF distribution
    so we don't need to correct the scaling factor here. When using
    LoadEventNexus, we have to scale by (t_frame-t_low_cut-t_high_cut)/t_frame.

    :param data: events workspace for data
    :param dark: events workspace for dark current
    :param logname: Log entry to calculate for duration. If None, duration will
    be tried looking sequentially into log entries duration', 'proton_charge',
    and 'timer.
    :return: Matrix workspace
    """
    ratio = duration_ratio(data.run(), dark.run(), logname=logname)
    return data - ratio * dark


class EQSANSDarkCurrentSubtract(PythonAlgorithm):

    def __init__(self):
        super(self.__class__, self).__init__()

    def category(self):
        return 'Workflow\\SANS\\EQSANS'

    def name(self):
        return 'EQSANSDarkCurrentSubtract'

    def summary(self):
        return "Perform EQSANS dark current subtraction"

    def version(self):
        return 1

    def PyInit(self):
        helpd = dict(Data='Events Workspace containing data from run(s)',
                     DarkCurrent='Events Workspace of the dark current',
                     LogName='Log entry to retrieve the duration. ' +
                     'If empty, the algorithm will do a sequential search ' +
                     'for entries "duration", "proton_charge", and "timer".',
                     OutputWorkspace='Subtracted events Workspace')

        def mwp(key, direction=Direction.Input):
            self.declareProperty(WorkspaceProperty(key, helpd[key], direction))
        mwp('Data')
        mwp('DarkCurrent')
        self.declareProperty('LogName', '', helpd['LogName'])
        mwp('OutputWorkspace', direction=Direction.Output)

    def PyExec(self):
        # Process input workspaces
        data = self.getProperty('Data').value
        dark = self.getProperty('DarkCurrent').value
        if data.id() != "EventWorkspace" or dark.id() != "EventWorkspace":
            raise AttributeError('Input workspace must be Event Workspaces')

        # Process input log entry name
        logname = self.getProperty('LogName').value
        if logname == '':
            logname = None

        subtracted = subtract_scaled_dark(data, dark, logname=logname)
        self.setProperty("OutputWorkspace", subtracted)

# Register algorithm with Mantid.
AlgorithmFactory.subscribe(EQSANSDarkCurrentSubtract)
