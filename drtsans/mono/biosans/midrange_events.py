# package imports
from drtsans.instruments import fetch_idf

# third party imports
from mantid.api import mtd
from mantid.simpleapi import LoadInstrument

# standard library imports
import tempfile


def has_midrange_detector(input_workpace):
    workspace = mtd[str(input_workpace)]
    return workspace.getInstrument().getComponentByName("midrange_detector") is not None


def update_idf(input_workspace):
    r"""
    Download the latest BIOSANS IDF from the Mantid GitHub repository and load it into a workspace.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to load the instrument into.

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with the instrument loaded.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        idf = fetch_idf("BIOSANS_Definition.xml", output_directory=temp_dir)
        LoadInstrument(Workspace=input_workspace, Filename=idf, RewriteSpectraMap=True)
    return mtd[str(input_workspace)]
