from drtsans.instruments import extract_run_number, instrument_enum_name, InstrumentEnumName
from drtsans.path import exists as path_exists
from drtsans.samplelogs import SampleLogs
from drtsans.settings import amend_config
# https://docs.mantidproject.org/nightly/api/python/mantid/api/AnalysisDataServiceImpl.html
from mantid.simpleapi import mtd
# https://docs.mantidproject.org/nightly/algorithms/LoadEventNexus-v1.html
from mantid.simpleapi import LoadEventNexus

__all__ = ['load_events']


def load_events(run, data_dir=None, output_workspace=None, overwrite_instrument=True, output_suffix='', **kwargs):
    r"""
    Load an event Nexus file produced by the HFIR instruments at ORNL.

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

    if mtd.doesExist(output_workspace):
        # if it exists skip loading
        return mtd[output_workspace]
    else:
        # load the data into the appropriate workspace
        with amend_config({'default.instrument': str(instrument_unique_name)}, data_dir=data_dir):
            # not loading the instrument xml from the nexus file will use the correct one that is inside mantid
            kwargs['LoadNexusInstrumentXML'] = not overwrite_instrument
            if 'LoadMonitors' not in kwargs:
                # load monitors for biosans or gpsans. eqsans does not.
                kwargs['LoadMonitors'] = is_mono
            LoadEventNexus(Filename=filename, OutputWorkspace=output_workspace, **kwargs)

    # insert monitor counts for monochromatic
    if is_mono:
        monitor_workspace = mtd[output_workspace + '_monitors']
        SampleLogs(output_workspace).insert('monitor', monitor_workspace.getNumberEvents())
        monitor_workspace.delete()

    return mtd[output_workspace]
