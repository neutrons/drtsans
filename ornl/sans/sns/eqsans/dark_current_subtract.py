from __future__ import (absolute_import, division, print_function)

from mantid.api import (PythonAlgorithm, AlgorithmFactory,
                        MatrixWorkspaceProperty)
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
        helpd = dict(Data='Matrix Workspace containing data from run(s)',
                     DarkCurrent='Matrix Workspace of the dark current',
                     OutputWorkspace='Subtracted matrix Workspace')

        def mwp(key, direction=Direction.Input):
            self.declareProperty(MatrixWorkspaceProperty(key, helpd[key],
                                                         direction))
        mwp('Data')
        mwp('DarkCurrent')
        mwp('OutputWorkspace', direction=Direction.Output)

    def PyExec(self):
        data = self.getProperty('Data').value
        dark = self.getProperty('DarkCurrent').value
        ratio = duration_ratio(data.run(), dark.run())
        subtracted = Scale(data, Factor=ratio, Operation='Multiply')
        self.setProperty("OutputWorkspace", subtracted)

# Register algorithm with Mantid.
AlgorithmFactory.subscribe(EQSANSDarkCurrentSubtract)
