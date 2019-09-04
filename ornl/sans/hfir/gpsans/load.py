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


def load_histogram(filename, wavelength=None, wavelength_spread=None, sample_to_detector_distance=None, unit='mm',
                   output_workspace=None):
    r"""

    Parameters
    ----------
    filename : string
        XML data file
    wavelength : float, optional
        The wavelength value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file,
        by default None
    wavelength_spread : float, optional
        wavelength spread value to use when loading the data file (Angstrom).
        This value will be used instead of the value found in the data file,
        by default None
    sample_to_detector_distance : float, optional
        Sample to detector distance to use (overrides meta data), in units of keyword `units`.
    unit: str
        Units for `sample_to_detector_distance`. Default is in mili-meters
    output_workspace : string
        The name of the output workspace. If none the name of the input file
        is used, by default None

    Returns
    -------
    MatrixWorkspace
        A reference to the workspace created.
    """
    unit_to_mm = dict(m=1.e3, mm=1.)
    if sample_to_detector_distance is not None:
        sample_to_detector_distance *= unit_to_mm[unit]
    if output_workspace is None:
        output_workspace = os.path.basename(filename).split('.')[0]
    LoadHFIRSANS(Filename=filename,
                 Wavelength=wavelength,
                 WavelengthSpread=wavelength_spread,
                 SampleDetectorDistance=sample_to_detector_distance,
                 OutputWorkspace=output_workspace)
    return mtd[output_workspace]
