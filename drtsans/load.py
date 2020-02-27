from drtsans.geometry import translate_detector_by_z, translate_sample_by_z
from drtsans.instruments import extract_run_number, instrument_enum_name, InstrumentEnumName
from drtsans.path import abspath
from drtsans.path import exists as path_exists
from drtsans.samplelogs import SampleLogs
from drtsans.settings import amend_config
import h5py
# https://docs.mantidproject.org/nightly/api/python/mantid/api/AnalysisDataServiceImpl.html
from mantid.simpleapi import mtd
# https://docs.mantidproject.org/nightly/algorithms/LoadEventNexus-v1.html
from mantid.simpleapi import LoadEventNexus

__all__ = ['load_events']


def __monitor_counts(filename, monitor_name='monitor1'):
    '''Get the total number of counts in a single monitor

    Parameters
    ----------
    filename: str
        Absolute path to file to be read
    monitor_name: str
        Name of the monitor to determine the total counts of
    '''
    counts = 0  # default value is zero
    with h5py.File(filename, 'r') as handle:
        if monitor_name not in handle['entry']:
            raise RuntimeError('File "{}" does not contain /entry/{}'.format(filename, monitor_name))
        # open the monitor group
        nxmonitor = handle['entry'][monitor_name]

        # get the number of counts from the total counts array or the monitor array
        if 'total_counts' in nxmonitor:
            counts = nxmonitor['total_counts'][0]
        else:
            counts = nxmonitor['event_time_offset'].shape[0]
    return int(counts)


def load_events(run, data_dir=None, output_workspace=None, overwrite_instrument=True, output_suffix='',
                detector_offset=0., sample_offset=0.,
                reuse_workspace=False, **kwargs):
    r"""
    Load an event Nexus file produced by the instruments at ORNL.

    Parameters
    ----------
    run: str, ~mantid.api.IEventWorkspace
        Examples: ``CG3_55555``, ``CG355555`` or file path.
    output_workspace: str
        If not specified it will be ``BIOSANS_55555`` determined from the supplied value of ``run``.
    data_dir: str, list
        Additional data search directories
    overwrite_instrument: bool, str
        If not :py:obj:`False`, ignore the instrument embedeed in the Nexus file. If :py:obj:`True`, use the
        latest instrument definition file (IDF) available in Mantid. If ``str``, then it should be the filepath to the
        desired IDF.
    output_suffix: str
        If the ``output_workspace`` is not specified, this is appended to the automatically generated
        output workspace name.
    detector_offset: float
        Additional translation of the detector along the Z-axis, in mm. Positive
        moves the detector downstream.
    sample_offset: float
        Additional translation of the sample, in mm. The sample flange remains
        at the origin of coordinates. Positive moves the sample downstream.
    reuse_workspace: bool
        When true, return the ``output_workspace`` if it already exists
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    instrument_unique_name = instrument_enum_name(run)  # determine which SANS instrument
    run_number = extract_run_number(run) if isinstance(run, str) else ''
    filename = run if path_exists(run) else '{}{}'.format(instrument_unique_name, run_number)

    # create default name for output workspace
    if (output_workspace is None) or (not output_workspace) or (output_workspace == 'None'):
        output_workspace = '{}_{}{}'.format(instrument_unique_name, run_number, output_suffix)

    # determine if this is a monochromatic measurement
    is_mono = (instrument_unique_name == InstrumentEnumName.BIOSANS) or \
              (instrument_unique_name == InstrumentEnumName.GPSANS)

    if reuse_workspace and mtd.doesExist(output_workspace):
        # if it exists skip loading
        return mtd[output_workspace]
    else:
        # load the data into the appropriate workspace
        with amend_config({'default.instrument': str(instrument_unique_name)}, data_dir=data_dir):
            # not loading the instrument xml from the nexus file will use the correct one that is inside mantid
            kwargs['LoadNexusInstrumentXML'] = not overwrite_instrument
            LoadEventNexus(Filename=filename, OutputWorkspace=output_workspace, **kwargs)

    # insert monitor counts for monochromatic instruments
    if is_mono:
        # determine the fully qualified file path
        if 'Filename' in mtd[output_workspace].run():
            # from the existing workspace
            filename = str(mtd[output_workspace].run()['Filename'].value)
        else:
            # use archive search
            filename = str(abspath(filename))

        # create new log with the monitor counts
        SampleLogs(output_workspace).insert('monitor', __monitor_counts(filename))

    # move instrument components - sample position must happen first
    translate_sample_by_z(output_workspace, 1e-3 * float(sample_offset))
    translate_detector_by_z(output_workspace, None)  # search logs and translate if necessary
    translate_detector_by_z(output_workspace, 1e-3 * float(detector_offset))

    return mtd[output_workspace]
