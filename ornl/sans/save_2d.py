from __future__ import print_function

from mantid.simpleapi import SaveNISTDAT, SaveNexus

def save_nist_dat(wksp, filename):
    """Save I(Qx, Qy) data to a text file compatible with NIST and DANSE readers

    Parameters
    ----------
    wksp : Workspace2D
        Workspace to be saved
    filename : string
        Filename of the output text file. Allowed extensions: ['.dat']
    """
    SaveNISTDAT(InputWorkspace=wksp, Filename=filename)


def save_nexus(wksp, title, filename):
    """Write the given Mantid workspace to a NeXus file.

    Parameters
    ----------
    wksp : Workspace2D
        Name of the workspace to be saved
    title : string
        Title to describe the saved worksapce
    filename : string
        The bame of the NeXus file to write, as a full or relative path. Allowed extensions: ['.nxs','.nx5','.xml']
    """
    SaveNexus(InputWorkspace=wksp, Title=title, Filename=filename)

