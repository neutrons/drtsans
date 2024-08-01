# package imports
from drtsans.instruments import fetch_idf

# third party imports
from mantid.api import mtd
from mantid.simpleapi import LoadInstrument

# standard library imports
import tempfile


def update_idf(input_workspace, **load_instrument_kwargs):
    r"""
    Download the latest BIOSANS IDF from the Mantid GitHub repository and load it into a workspace.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to load the instrument into.

    load_instrument_kwargs : dict
        Keyword arguments to pass to Mantid algorithm LoadInstrument.

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with the instrument loaded.
    """
    kwargs = dict(RewriteSpectraMap=False)  # default options to Mantid algorithm LoadInstrument
    kwargs.update(load_instrument_kwargs)

    with tempfile.TemporaryDirectory() as temp_dir:
        idf = fetch_idf("BIOSANS_Definition.xml", output_directory=temp_dir)
        LoadInstrument(Workspace=input_workspace, Filename=idf, **kwargs)
    return mtd[str(input_workspace)]
