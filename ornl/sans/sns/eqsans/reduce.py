from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import Load

from ornl.settings import (amend_config, optional_output_workspace,
                           unique_workspace_dundername as uwd)
from ornl.sans import geometry
from ornl.sans.sns.eqsans import geometry as e_geometry, correct_frame


@optional_output_workspace
def load_w(run, low_tof_clip=0, high_tof_clip=0, dw=0.1):
    r"""
    Load a run, correct the TOF frame, and convert to wavelength

    Parameters
    ----------
    run: str
        Run number or filename. Passed onto Mantid's `Load` algorithm
    low_tof_clip: float
        Lower TOF clipping
    high_tof_clip: float
        Upper TOF clipping
    dw: float
        Wavelength bin width

    Returns
    -------
    MatrixWorkspace
    """
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        _ws = Load(run, OutputWorkspace=uwd())
        e_geometry.translate_detector_z(_ws)
        correct_frame.correct_detector_frame(_ws)
        sdd = geometry.source_detector_distance(_ws, units='m')
        args = (_ws, sdd, low_tof_clip, high_tof_clip)
        bands = correct_frame.transmitted_bands_clipped(*args,
                                                        interior_clip=True)
        _ws = correct_frame.convert_to_wavelength(_ws, bands, dw, events=False)
        correct_frame.log_tof_structure(_ws, low_tof_clip, high_tof_clip,
                                        interior_clip=True)
        return _ws
