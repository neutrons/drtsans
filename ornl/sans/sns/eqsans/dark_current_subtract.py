from __future__ import (absolute_import, division, print_function)

from mantid.api import (PythonAlgorithm, AlgorithmFactory, WorkspaceProperty)
from mantid.simpleapi import (SumSpectra, AppendSpectra, ExtractSpectra,
                              CompressEvents)
from mantid.kernel import Direction, logger, StringListValidator

from ornl.sans.mask_utils import masked_indexes


def compute_log_ratio(run1, run2, l):
    """Compute ratio of data to dark for one log entry
    :param run1: first run object containing logs
    :param run2: second run object containing logs
    :param l: entry log
    :return: ratio
    :except: RuntimeError when log is not present
    """
    dt = run1.getProperty(l)
    dk = run2.getProperty(l)
    try:
        dt = dt.getStatistics().duration
        dk = dk.getStatistics().duration
    except AttributeError:
        dt = dt.value
        dk = dk.value
    return dt / dk


def duration_ratio(data, dark, log_name=None):
    """Compute the ratio of data to dark-current durations
    :param data: run object for data
    :param darkcurrent: run object for dark current
    :param log_name: entry log containing the duration. If None, duration will
    be tried looking sequentially into log entries duration', 'proton_charge',
    and 'timer.
    :return: duration ratio, or 1.0 if ratio cannot be computed
    """
    not_found = 'Logs could not be found, duration ratio set to 1.0'
    if log_name is not None:
        try:
            return compute_log_ratio(data, dark, log_name)
        except RuntimeError:
            logger.error(not_found)
    else:
        for l in ('duration', 'proton_charge', 'timer'):
            try:
                return compute_log_ratio(data, dark, l)
            except RuntimeError:
                continue
    return 1.0


def subtract_pixelcount_dark(data, dark, log_name=None):
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
    :param log_name: Log entry to calculate for duration. If None, duration will
    be tried looking sequentially into log entries 'duration', 'proton_charge',
    and 'timer'
    :return: events workspace
    """
    ratio = duration_ratio(data.run(), dark.run(), log_name=log_name)
    return data - ratio * dark


def subtract_isotropic_dark(data, dark, log_name=None):
    """Integrate dark counts, rescale, and subtract from data

    All dark events in unmasked pixels are summed up, rescaled by the
    duration_ration and the number of unmasked pixels. This provides an
    effective dark event list that can be applied to each data pixel.

    The scaling factor should account for the TOF cuts on each side of a frame
    The EQSANSLoad algorithm cuts the beginning and end of the TOF distribution
    so we don't need to correct the scaling factor here. When using
    LoadEventNexus, we have to scale by (t_frame-t_low_cut-t_high_cut)/t_frame.

    :param data: events workspace for data
    :param dark: events workspace for dark current
    :param log_name: Log entry to calculate for duration. If None, duration will
    be tried looking sequentially into log entries 'duration', 'proton_charge',
    and 'timer'
    :return: events workspace
    """
    # Collect all dark events listed in the unmasked pixels into a single
    # event list, then rescale by number of unmasked pixels to yield an
    # effective list of dark events that can be applied to a single pixel
    dark_unmasked_indexes = masked_indexes(dark, invert=True).tolist()
    dark_summed = SumSpectra(dark,
                             ListOfWorkspaceIndices=dark_unmasked_indexes)
    dark_summed /= len(dark_unmasked_indexes)
    # Rescale list of dark events
    dark_summed *= duration_ratio(data.run(), dark.run(), log_name=log_name)
    # Compress, to avoid an unmanageable number of dark events
    dark_events_per_pixel = 1000  # we aim to this many dark events per pixel
    #   Add entries for other units, like wavelength, if needed
    initial_tolerances = {'Time-of-flight': 1.0,  # TOF in microseconds
                          }
    units = dark_summed.getAxis(0).getUnit().name()
    tolerance = initial_tolerances[units]
    while dark_summed.getNumberEvents() > dark_events_per_pixel:
        dark_summed = CompressEvents(dark_summed, Tolerance=tolerance)
        tolerance *= 2.0
    # Repeat the dark event lists for every histogram in data workspace
    dark_summed.getSpectrum(0).clearDetectorIDs()
    n_data_histograms = data.getNumberHistograms()
    while dark_summed.getNumberHistograms() < n_data_histograms:
        dark_summed = AppendSpectra(dark_summed, dark_summed,
                                    ValidateInputs=False, MergeLogs=False)
    dark_replicated = ExtractSpectra(dark_summed,
                                     EndWorkspaceIndex=n_data_histograms-1)
    return data - dark_replicated


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
                     Method='Subtraction algorithm',
                     OutputWorkspace='Subtracted events Workspace')

        def mwp(key, direction=Direction.Input):
            self.declareProperty(WorkspaceProperty(key, helpd[key], direction))
        mwp('Data')
        mwp('DarkCurrent')
        self.declareProperty('LogName', '', helpd['LogName'])
        method_validator = StringListValidator(['PixelCount', 'Isotropic'])
        self.declareProperty('Method', defaultValue='PixelCount',
                             validator=method_validator,
                             direction=Direction.Input, doc=helpd['Method'])
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

        method = self.getProperty('Method').value
        method_func = dict(PixelCount=subtract_pixelcount_dark,
                           Isotropic=subtract_isotropic_dark)
        subtracted = method_func[method](data, dark, logname=logname)

        self.setProperty("OutputWorkspace", subtracted)

# Register algorithm with Mantid.
AlgorithmFactory.subscribe(EQSANSDarkCurrentSubtract)
