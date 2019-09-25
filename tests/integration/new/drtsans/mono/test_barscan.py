import pytest
import numpy as np
from drtsans.mono.barscan import find_edges, fit_positions

r"""Finding the edges of the barscan in a single tube,
then calculate the position and width of the pixels

based on BarScanShadow_test_KCL_SVP.xlsx and BarScanFitTest_KCL.xlsx
Testing appendix 2.1 in the master document

devs - Andrei Savici <saviciat@ornl.gov>,
       Jose Borreguero <borreguerojm@ornl.gov>
SME  - Ken Littrell <littrellkc@ornl.gov>
"""


def test_find_edges():
    intensities = np.array([2, 2, 38, 38, 38, 34, 38, 41, 35, 3, 4, 3, 3, 4,
                            30, 30, 37, 33, 31, 39, 42, 42, 2, 2, 2])
    edges = find_edges(intensities)
    assert (edges.bottom_pixel, edges.top_pixel) == (2, 21)
    assert (edges.bottom_shadow_pixel, edges.above_shadow_pixel) == (9, 14)
    with pytest.raises(RuntimeError, match='Faulty tube found'):
        find_edges(intensities, min_illuminated_length=len(intensities))


def test_no_shadow():
    intensities = [1]*25
    with pytest.raises(IndexError, match='Could not find bottom shadow edge'):
        find_edges(intensities)


def test_no_bottom_tube_edge():
    intensities = np.ones(25)
    intensities[14] = 3000  # all pixels but one are below the threshold
    with pytest.raises(IndexError, match='Could not find bottom tube edge'):
        find_edges(intensities)


def test_fit_positions():
    edge_pixel = np.arange(25) + 1
    pos = [-120, -110, -98, -87, -74, -62, -49, -36, -23, -11,
           1, 12, 23, 32, 41, 49, 57, 64, 71, 78, 86, 93, 102, 110, 120]
    new_positions, new_heights = fit_positions(edge_pixel, pos, tube_pixels=26)
    expected_positions = [-119.659,
                          -109.907,
                          -98.892,
                          -86.953,
                          -74.396,
                          -61.496,
                          -48.494,
                          -35.599,
                          -22.989,
                          -10.808,
                          0.8303,
                          11.846,
                          22.189,
                          31.846,
                          40.831,
                          49.193,
                          57.013,
                          64.403,
                          71.508,
                          78.506,
                          85.606,
                          93.049,
                          101.109,
                          110.094,
                          120.340]
    expected_heights = [9.001,
                        10.443,
                        11.530,
                        12.296,
                        12.771,
                        12.989,
                        12.981,
                        12.779,
                        12.417,
                        11.926,
                        11.337,
                        10.685,
                        10.000,
                        9.315,
                        8.663,
                        8.075,
                        7.583,
                        7.221,
                        7.019,
                        7.011,
                        7.228,
                        7.704,
                        8.469,
                        9.556,
                        10.998]
    # fit_positions calculates also the expected position for pixel 0, not in the table
    assert new_positions[1:] == pytest.approx(expected_positions, abs=1e-2)
    assert new_heights[1:] == pytest.approx(expected_heights, abs=1e-2)


if __name__ == '__main__':
    pytest.main()
