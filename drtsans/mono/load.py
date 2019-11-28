import os

# https://docs.mantidproject.org/nightly/algorithms/LoadHFIRSANS-v1.html
from mantid.simpleapi import LoadHFIRSANS, LoadEventNexus, CloneWorkspace, LoadInstrument, HFIRSANS2Wavelength
from mantid.api import mtd

from drtsans.instruments import extract_run_number, instrument_enum_name
from drtsans.settings import amend_config
from drtsans.samplelogs import SampleLogs
from drtsans.process_uncertainties import set_init_uncertainties

__all__ = ['load_events', 'load_histogram', 'transform_to_wavelength', 'load_mono']


def load_events(run, output_workspace=None, data_dir=None, overwrite_instrument=False, **kwargs):
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
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    instrument_unique_name = str(instrument_enum_name(run))  # elucidate which SANS instrument
    run_number = str(extract_run_number(run)) if isinstance(run, str) else ''

    if output_workspace is None:
        output_workspace = instrument_unique_name + '_' + run_number

    if isinstance(run, str):
        with amend_config({'default.instrument': instrument_unique_name}, data_dir=data_dir):
            if overwrite_instrument is not False:
                kwargs['LoadNexusInstrumentXML'] = False
            LoadEventNexus(Filename=instrument_unique_name + run_number, OutputWorkspace=output_workspace,
                           LoadMonitors=True, **kwargs)
    else:
        CloneWorkspace(run, OutputWorkspace=output_workspace)

    if overwrite_instrument is not False:
        optional_arguments = dict(InstrumentName=instrument_unique_name) if overwrite_instrument is True else dict(
            Filename=overwrite_instrument)
        LoadInstrument(Workspace=output_workspace, RewriteSpectraMap=False, **optional_arguments)

    # Insert monitor counts as log entry
    monitor_workspace = mtd[output_workspace + '_monitors']
    SampleLogs(output_workspace).insert('monitor', monitor_workspace.getNumberEvents())
    monitor_workspace.delete()

    return mtd[output_workspace]


def load_histogram(filename, output_workspace=None, wavelength=None, wavelength_spread=None, sample_det_cent=None):
    """Loads a SANS data file produce by the HFIR instruments at ORNL.
    The instrument geometry is also loaded. The center of the detector is
    placed at (0, 0, :ref:`sample_det_cent <devdocs-standardnames>` )

    Parameters
    ----------
    filename : str
        The name of the input xml file to load
    output_workspace : str, optional
        The optional name of the output workspace. If :py:obj:`None` is the filename stripped of the extension.
    wavelength : float
        The wavelength value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file.
    wavelength_spread : float
        wavelength spread value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file.
    sample_det_cent : float
        Sample to detector distance to use (overrides meta data) in mm

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        A reference for the workspace created.
    """

    if output_workspace is None:
        output_workspace = os.path.basename(filename).split('.')[0]

    ws = LoadHFIRSANS(Filename=filename, Wavelength=wavelength, WavelengthSpread=wavelength_spread,
                      SampleDetectorDistance=sample_det_cent, OutputWorkspace=output_workspace)
    return ws


def transform_to_wavelength(input_workspace, output_workspace=None):
    r"""
    Transforms the event files with fake time of flight from the SANS instruments at HFIR into histograms
    in wavelength.

    **Mantid Algorithms used:**
    :ref:`Divide <algm-HFIRSANS2Wavelength-v1>`,

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace
        Events workspace in time-of-flight.
    output_workspace: str
        Name of the output workspace. If :py:obj:`None`, the name of the input_workspace will be
        used, thus overwritting the input workspace.
    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    HFIRSANS2Wavelength(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)

    # Set initial uncertainties
    input_workspace = set_init_uncertainties(output_workspace)

    return mtd[output_workspace]


def load_mono(filename, **kwargs):
    r"""
    Loads a SANS data file produce by the HFIR instruments at ORNL.

    Parameters
    ----------

    filename: int, str
        Examples: ``55555`` or ``CG3_55555`` or file path.
    kwargs:
        keyword arguments for load_events or load_histogram.
    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    try:
        return load_events(filename, **kwargs)
    except Exception:
        return load_histogram(filename, **kwargs)
