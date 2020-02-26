"""
    Test EASANS sensitivities preparation algorithm
"""
import pytest
import numpy as np
from drtsans.mono.gpsans import prepare_data, find_beam_center


def test_gpsans_find_beam_center():
    """Integration test on algorithm to find beam center for GPSANS

    Returns
    -------

    """
    # Load data
    beam_center_ws = prepare_data(data='CG2_8148',
                                  btp={'Pixel': '1-8,249-256'},
                                  detector_offset=0,
                                  sample_offset=0,
                                  center_x=0,
                                  center_y=0,
                                  flux_method='monitor',
                                  solid_angle=False,
                                  output_workspace='BC_8148')

    # Find beam center
    beam_center = find_beam_center(beam_center_ws)

    # Get detector center
    instrument = beam_center_ws.getInstrument()
    det = instrument.getComponentByName('detector1')
    det_center = det.getPos()

    # Calculate shift:
    center_x, center_y = beam_center
    beam_center_shift = np.sqrt((center_x - det_center[0])**2 + (center_y - det_center[1])**2)

    assert beam_center_shift == pytest.approx(0.400, abs=0.007), 'Beam center shift {} to {} is beyond' \
                                                                 '0.4 +/- 7E-3'.format(beam_center,
                                                                                       det_center)


if __name__ == '__main__':
    pytest.main([__file__])
