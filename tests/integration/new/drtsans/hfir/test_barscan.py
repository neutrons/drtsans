import pytest
import numpy as np
from ornl.sans.hfir.barscan import find_edges

r"""Finding the edges of the barscan in a single tube
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


if __name__ == '__main__':
    pytest.main()
