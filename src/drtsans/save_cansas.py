"""
Methods for CANSAS format handling
"""

from mantid.kernel import logger
from mantid.simpleapi import SaveCanSAS1D, SaveNXcanSAS


def save_cansas_nx(*args, **kwargs):
    """Writes a MatrixWorkspace to a file in the NXcanSAS format

    Parameters
    ----------
    InputWorkspace : MatrixWorkspace
        The input workspace, which must be in units of Q
    filename : string
        The name of the .h5 file to save. Allowed extensions: [‘.h5’]

    See more: https://docs.mantidproject.org/v6.1.0/algorithms/SaveNXcanSAS-v1.html
    """

    try:
        SaveNXcanSAS(*args, **kwargs)
    except ValueError as e:
        if "The workspace must have common bin boundaries for all histograms" in str(e):
            logger.warning(str(e))
        else:
            raise e


def save_cansas_xml_1D(wksp, title, filename):
    """Save the I(q) workspace in SaveCanSAS (XML) format

    Parameters
    ----------
    wksp : ~mantid.api.MatrixWorkspace
        Workspace containing only one spectrum (the I(q) curve)
    title : string
        Text to append to Process section
    filename : string
        The output filename
    """
    SaveCanSAS1D(InputWorkspace=wksp, Process=title, Filename=filename)
