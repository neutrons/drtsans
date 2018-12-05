from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import Load

from ornl.settings import amend_config, uwn
from ornl.sans import geometry
from ornl.sans.sns.eqsans import geometry as e_geometry, correct_frame


def load_w(run, out_ws, ltc=0, htc=0, dw=0.1):
    r"""
    Load a run, correct the TOF frame, and convert to wavelength

    Parameters
    ----------
    run: str
        Run number or filename. Passed onto Mantid's `Load` algorithm
    out_ws: str
        Name of the output workspace
    ltc: float
        Lower TOF clipping
    htc: float
        Upper TOF clipping
    dw: float
        Wavelength bin width

    Returns
    -------
    MatrixWorkspace
    """
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        _ws = Load(run, OutputWorkspace=uwn())
        try:
            e_geometry.translate_detector_z(_ws)
            correct_frame.correct_detector_frame(_ws)
            sdd = geometry.source_detector_distance(_ws, units='m')
            bands = correct_frame.transmitted_bands_clipped(_ws, sdd, ltc, htc,
                                                            interior_clip=True)
            wl = correct_frame.convert_to_wavelength(_ws, bands, dw, out_ws,
                                                     events=False)
            correct_frame.log_tof_structure(wl, ltc, htc, interior_clip=True)
        finally:
            _ws.delete()
        return wl
