import os

# https://docs.mantidproject.org/nightly/algorithms/LoadHFIRSANS-v1.html
from mantid.simpleapi import LoadHFIRSANS, LoadEventNexus, CloneWorkspace
from mantid.api import mtd

from drtsans.instruments import instrument_enum_name
from drtsans.samplelogs import SampleLogs

__all__ = ['load_events', 'load_histogram']


def load_events(run, output_workspace=None, data_dir=None, **kwargs):
    r"""
    Load an event Nexus file produced by the HFIR instruments at ORNL.

    Parameters
    ----------
    run: int, str
        Examples: ``55555`` or ``CG3_55555`` or file path.
    data_dir: str, list
        Additional data search directories
    output_workspace: str
        If not specified it will be ``BIOSANS_55555`` determined from the supplied value of ``run``.
    kwargs: dict
        Additional positional arguments for :ref:`LoadEventNexus <algm-LoadEventNexus-v1>`.

    Returns
    -------
    ~mantid.api.IEventWorkspace
        Reference to the events workspace
    """
    instrument_name = str(instrument_enum_name(run))  # elucidate which SANS instrument
    if output_workspace is None:
        if isinstance(run, str):
            output_workspace = os.path.split(run)[-1]
            output_workspace = instrument_name + '_' + output_workspace.split('_')[1]
            output_workspace = output_workspace.split('.')[0]
        else:
            output_workspace = instrument_name + '_' + str(run)

    if isinstance(run, int) or isinstance(run, str):
        with amend_config({'default.instrument': instrument_name}, data_dir=data_dir):
            LoadEventNexus(Filename=str(run), OutputWorkspace=output_workspace, LoadMonitors=True, **kwargs)
    else:
        CloneWorkspace(run, OutputWorkspace=output_workspace)

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
