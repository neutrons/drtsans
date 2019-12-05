import numpy as np
import pytest

r""" Hyperlinks to mantid algorithms
LoadEmptyInstrument <https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html>
"""

r"""
Hyperlinks to drtsans functions
ElementComponentInfo, PixelInfo, TubeInfo, TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.settings import unique_workspace_dundername
from drtsans.tubecollection import TubeCollection
from drtsans.mono.biosans import load_events


@pytest.fixture(scope='module')
def tube_collection(reference_dir):
    workspace = load_events('CG3_961.nxs.h5',
                            output_workspace=unique_workspace_dundername(),
                            data_dir=reference_dir.new.biosans,
                            overwrite_instrument=True)
    collection = TubeCollection(workspace, 'detector1')
    return collection


@pytest.mark.usefixtures('tube_collection')
class TestTubeCollection(object):

    def test_tubes(self, tube_collection):
        assert len(tube_collection.tubes) == 192

    def test_getitem(self, tube_collection):
        assert tube_collection.tubes[0][0].position == pytest.approx([0.53, -0.52, 3.05], abs=0.01)  # first pixel
        assert tube_collection.tubes[-1][-1].position == pytest.approx([-0.54, 0.52, 3.06], abs=0.01)  # last pixel

    def test_sorted(self, tube_collection):
        # Sort by decreasing tube position along the X-axis
        sorted_tubes = tube_collection.sorted(view='decreasing X')
        x_coords = [tube[0].position[0] for tube in sorted_tubes]  # X coord for the first pixel of each tube
        assert np.all(x_coords[1:] < x_coords[:-1])  # x_coords strictly decreasing
        # Sort by increasing spectrum index
        sorted_tubes = tube_collection.sorted(view='spectrum index')
        spectrum_indexes = [tube[0].spectrum_index for tube in sorted_tubes]
        assert np.all(spectrum_indexes[1:] > spectrum_indexes[:-1])


if __name__ == '__main__':
    pytest.main([__file__])
