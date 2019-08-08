import os

from mantid.simpleapi import LoadHFIRSANS
from ornl.sans import solid_angle_correction


def apply_solid_angle_correction(input_workspace):
    """
        Apply solid angle correction
    """
    return solid_angle_correction(input_workspace, detector_type='VerticalTube')


def load(filename, output_workspace=None, wavelength=None, wavelength_spread=None,
         sample_to_detector_distance=None):
    """Loads a SANS data file produce by the HFIR instruments at ORNL.
    The instrument geometry is also loaded. The center of the detector is 
    placed at (0,0,D), where D is the sample-to-detector distance.

    Parameters
    ----------
    filename : string
        The name of the input xml file to load
        output_workspace : string, optional
        The name of the Output workspace. If none is the filename stripped
        of the extension, by default None
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
    """

    if output_workspace is None:
        output_workspace = os.path.basename(filename).split('.')[0]

    ws = LoadHFIRSANS(filename, output_workspace, Wavelength=wavelength,
                      WavelengthSpread=wavelength_spread,
                      SampleDetectorDistance=sample_to_detector_distance)
    return ws
