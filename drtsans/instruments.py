# local imports
from drtsans.settings import unique_workspace_dundername

# third party imports
from mantid.kernel import ConfigService
from mantid.api import mtd, MatrixWorkspace
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import LoadInstrument, LoadEmptyInstrument, MergeRuns, RemoveSpectra, RenameWorkspace

# standard imports
import enum
import os
import subprocess
from typing import Optional, Union


__all__ = [
    "InstrumentEnumName",
    "instrument_enum_name",
    "instrument_standard_name",
    "is_time_of_flight",
]

INSTRUMENT_LABELS = ["CG3", "BIOSANS", "EQ-SANS", "EQSANS", "CG2", "GPSANS"]


@enum.unique
class InstrumentEnumName(enum.Enum):
    @staticmethod
    def names():
        r"""Standard names for all instruments, in alphabetical order"""
        names_all = list(map(str, InstrumentEnumName))
        names_all.remove("UNDEFINED")
        return sorted(names_all)

    r"""Unique names labelling each instrument"""
    UNDEFINED = None  # usually the dummy instrument used for testing
    BIOSANS = ConfigService.getFacility("HFIR").instrument("BIOSANS")
    EQSANS = ConfigService.getFacility("SNS").instrument("EQSANS")
    GPSANS = ConfigService.getFacility("HFIR").instrument("GPSANS")

    def __str__(self):
        return self.name


def instrument_enum_name(input_query):
    r"""
    Resolve the instrument name as a unique enumeration.

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    InstrumentEnumName
        The name of the instrument as one of the InstrumentName enumerations
    """
    string_to_enum = {
        "CG3": InstrumentEnumName.BIOSANS,
        "BIOSANS": InstrumentEnumName.BIOSANS,
        "EQ-SANS": InstrumentEnumName.EQSANS,
        "EQSANS": InstrumentEnumName.EQSANS,
        "CG2": InstrumentEnumName.GPSANS,
        "GPSANS": InstrumentEnumName.GPSANS,
    }
    # convert to a string
    name = str(input_query)

    if name in mtd:  # convert mantid workspace into a instrument string
        name = mtd[str(name)].getInstrument().getName()
    else:  # see if `name` contains any of the instrument labels
        name = name.upper()
        for instrument_string_label in sorted(string_to_enum.keys()):
            if instrument_string_label in name:
                name = instrument_string_label
                break
    return string_to_enum.get(name.upper(), InstrumentEnumName.UNDEFINED)


def instrument_standard_name(input_query):
    r"""
    Resolve the standard instrument name.

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    str
        The name of the instrument as the string representation of one of the InstrumentName enumerations
    """
    return str(instrument_enum_name(input_query))


def instrument_standard_names():
    r"""Standard names for all instruments, in alphabetical order"""
    return InstrumentEnumName.names()


def instrument_filesystem_name(input_query):
    r"""
    Resolve the name of the instrument that is the subdirectory name under /SNS or /HFIR

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    str
    """
    filesystem_name = {"BIOSANS": "CG3", "EQSANS": "EQSANS", "GPSANS": "CG2"}
    return filesystem_name[instrument_standard_name(input_query)]


def instrument_label(input_query):
    r"""
    Resolve the instrument name.

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    str
    """
    # convert to a string
    name = str(input_query)

    if name in mtd:  # convert mantid workspace into a instrument string
        return mtd[str(name)].getInstrument().getName()
    else:  # see if `name` contains any of the instrument labels
        name = name.upper()
        for instrument_string_label in INSTRUMENT_LABELS:
            if instrument_string_label in name:
                return instrument_string_label
    raise RuntimeError('Instrument name can not be resolved from "{}"'.format(input_query))


def extract_run_number(input_query):
    r"""
    Extract the run number from string

    Example:
    input string '/HFIR/..../CG3_961.nxs.h5', 'CG3_961.nxs.h5', 'CG3961', and 'CG3_961' should all return run
    number 961

    Parameters
    ----------
    input_query: str

    Returns
    -------
    int
    """
    try:
        # see if `input_query` is an integer
        run_number = int(input_query)
    except ValueError:
        # name of the file without path
        run_number = os.path.basename(input_query)
        # everything up to the extension
        run_number = run_number.split(".")[0]
        # remove the instrument name
        for label in INSTRUMENT_LABELS:
            run_number = run_number.replace(label, "")
        # remove any remaining '_'
        if "_" in run_number:
            run_number = run_number.split("_")[1]
        # convert to an integer

    return int(run_number)


def is_time_of_flight(input_query):
    r"""
    Find if the instrument is a time-of-flight one

    Parameters
    ----------
    input_query: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace, InstrumentEnumName
        string representing a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    bool
    """
    return instrument_enum_name(input_query) is InstrumentEnumName.EQSANS  # we only have one, for the moment


def fetch_idf(idf_xml, output_directory=os.getcwd()):
    r"""
    Download an IDF from the Mantid GitHub repository to a temporary directory.

    Parameters
    ----------
    idf_xml : str
        The name of the IDF file to download.
    output_directory: str

    Returns
    -------
    str
        absolute path to the downloaded IDF file.
    """
    idf = os.path.join(output_directory, idf_xml)
    url = f"https://raw.githubusercontent.com/mantidproject/mantid/main/instrument/{idf_xml}"
    result = subprocess.run(f"curl -o {idf} {url}", shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Dowloading {idf_xml} failed with error: {result.stderr}")
    return idf


def empty_instrument_workspace(
    output_workspace: str,
    filename: Optional[str] = None,
    instrument_name: Optional[str] = None,
    event_workspace: Optional[bool] = False,
    monitors_have_spectra: Optional[bool] = False,
) -> Union[MatrixWorkspace, EventWorkspace]:
    r"""
    Create an emtpy workspace for one of the standard instruments. By default, monitors do not have associated spectra.

    Invokes ~mantid.simpleapi.LoadEmptyInstrument to create an empty instrument workspace.

    Parameters
    ----------
    output_workspace
        Name of the output workspace
    filename
        Path to the instrument filename. If not absolute path, mantid will search in the instrument directory.
    instrument_name
        Alternative of option ``filename``. Mantid will search for the latest instrument file in the instrument
         directory.
    event_workspace
        If True, create an event workspace, otherwise a histogram workspace.
    monitors_have_spectra
        If True, create a workspace with spectra for the monitors.

    Returns
    -------
    A handle to the empty workspace
    """
    workspace = LoadEmptyInstrument(
        OutputWorkspace=output_workspace,
        Filename=filename,
        InstrumentName=instrument_name,
        MakeEventWorkspace=event_workspace,
    )
    if monitors_have_spectra is False:
        # get the list of non-negative detector IDs (i.e. excluding monitors)
        detector_ids = workspace.detectorInfo().detectorIDs()
        monitor_count = detector_ids[detector_ids < 0].size  # monitors have negative IDs always
        detector_ids = detector_ids[detector_ids >= 0]  # these are the detector IDs for detector pixels
        # iterate over the spectra, assigning one detector ID to each spectrum
        for wi, detid in zip(range(workspace.getNumberHistograms()), detector_ids):
            spectrum = workspace.getSpectrum(wi)
            spectrum.setDetectorID(int(detid))
        # the number of spectra should be equal to the number of monitors plus the number detector pixels. Thus,
        # we have an excess of spectra in the amount of `monitor_count`. They need to be removed now
        to_remove = list(range(workspace.getNumberHistograms() - monitor_count, workspace.getNumberHistograms()))
        workspace = RemoveSpectra(
            Inputworkspace=output_workspace, OutputWorkspace=output_workspace, WorkspaceIndices=to_remove
        )
    return workspace


def copy_to_newest_instrument(
    input_workspace: Union[str, MatrixWorkspace, EventWorkspace],
    output_workspace: Optional[str] = None,
) -> Union[MatrixWorkspace, EventWorkspace]:
    r"""
    Copy the workspace intensities and/or events to the latest instrument file.

    Will also copy the logs and preserve the original locations of the main and
    wing detectors.

    Parameters
    ----------
    input_workspace
        Workspace containing an old instrument definition file (IDF).
    output_workspace
        Workspace containing the intensities and/or events of input_workspace but with the
         latest IDF. If ``None``, the input workspace is overwritten
    """
    origin = mtd[str(input_workspace)]
    if output_workspace is None:
        target_workspace = unique_workspace_dundername()  # temporary name
    else:
        target_workspace = output_workspace
    instrument_file = fetch_idf(f"{origin.getInstrument().getName()}_Definition.xml")
    target = empty_instrument_workspace(
        output_workspace=target_workspace,
        filename=instrument_file,
        event_workspace=isinstance(origin, EventWorkspace),
        monitors_have_spectra=(origin.getSpectrum(0).getDetectorIDs()[0] < 0),
    )
    # for algorithm MergeRuns to work, units of origin and target workspace must match
    origin_unit = origin.getAxis(0).getUnit().unitID()
    target.getAxis(0).setUnit(origin_unit)
    target.setYUnit(origin.YUnit())
    MergeRuns(
        InputWorkspaces=[target_workspace, input_workspace], OutputWorkspace=target_workspace  # order is necessary
    )
    # Move components to the positions they have in input_workspace by reading their positions
    # in the logs. This is implicitly done when invoking algorithm LoadInstrument.
    LoadInstrument(Workspace=target_workspace, Filename=instrument_file, RewriteSpectraMap=False)

    if output_workspace is None:  # overwrite the input workspace
        target = RenameWorkspace(InputWorkspace=target_workspace, OutputWorkspace=str(input_workspace))
    return target
