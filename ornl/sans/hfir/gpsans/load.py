import os
import re
import numpy as np
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from mantid.simpleapi import (mtd, LoadHFIRSANS, LoadInstrument,
                              MoveInstrumentComponent)
from ornl.path import exists

# Functions exposed to the general user (public) API
__all__ = ['load_histogram']


def permute_contents(filename):
    r"""
    Permute the rows in the XML data file to match the tube ordering of the
    new IDF.

    The XML data file assumes the following ordering of tube ID's within the
    following eightpack (seen from the top)
                          7o  5o  3o  1o
                            6o  4o  2o  0o
    The new IDF has a different ordering
                          7o  6o  5o  4o
                            3o  2o  1o  0o
    Thus a permutation of the indexes, corresponding a permutation of the
    rows in the XML data file, is necessary prior to applying the new IDF.

    Parameters
    ----------
    filename: str
        Original data filename
    Returns
    -------
    str
        contents of original with permuted data rows
    """
    # Permutation from new to old tube indexes within an eightpack
    perm = list(range(0, 8, 2)) + list(range(1, 8, 2))
    # Iterate over the 24 eightpacks
    old_ids = [i + 8 * e for e in range(24) for i in perm]

    # Extract the tube counts into rows and permute
    pattern = r'<Detector type="INT32\[192,256\]">\n([\s\S]+)\n</Detector>'
    with open(filename) as f:
        contents = f.read()
    old_counts = re.search(pattern, contents).group(1)
    rows = np.array(old_counts.split('\n'), dtype=object)[old_ids]

    return contents.replace(old_counts, '\n'.join(rows))


@contextmanager
def serve_data_file(filename, idf=None):
    r"""
    If needed, create a temporary data file suitable for the new IDF.

    Parameters
    ----------
    filename: str
        XML or nexus file containing the data
    idf: str
        File path to new instrument definition file overriding the default one.

    Returns
    -------
    str
    """
    data_file = filename
    if idf is not None and exists(idf) and filename[-3:] == 'xml':
        f = NamedTemporaryFile('w', delete=False)
        f.write(permute_contents(filename))
        data_file = f.name
    try:
        yield data_file
    finally:
        if data_file != filename:
            os.remove(data_file)


def load_histogram(filename, wavelength=None, wavelength_spread=None,
                   sample_to_detector_distance=None,
                   idf=None, output_workspace=None, ):
    r"""

    Parameters
    ----------
    filename : string
        XML or Nexus file
    wavelength : float, optional
        The wavelength value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file,
        by default None
    wavelength_spread : float, optional
        wavelength spread value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file,
        by default None
    sample_to_detector_distance : float, optional
        Sample to detector distance to use (overrides meta data) in mm,
        by default None
    idf: str
        File path to instrument definition file overriding the loaded one
    output_workspace : string
        The name of the output workspace. If none the name of the input file
        is used, by default None

    Returns
    -------
    MatrixWorkspace
        A reference to the workspace created.
    """

    if output_workspace is None:
        output_workspace = os.path.basename(filename).split('.')[0]
    with serve_data_file(filename, idf) as data_file:
        LoadHFIRSANS(Filename=data_file,
                     Wavelength=wavelength,
                     WavelengthSpread=wavelength_spread,
                     SampleDetectorDistance=sample_to_detector_distance,
                     OutputWorkspace=output_workspace)
    if idf is not None and exists(idf):
        LoadInstrument(Workspace=output_workspace, Filename=idf,
                       RewriteSpectraMap=True)
        z = mtd[output_workspace].getInstrument().getSample().getPos()[-1] +\
            sample_to_detector_distance / 1e3
        MoveInstrumentComponent(Workspace=output_workspace,
                                ComponentName='detector1',
                                Z=z, RelativePosition=False)
    return mtd[output_workspace]
