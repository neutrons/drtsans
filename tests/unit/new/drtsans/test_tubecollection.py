import numpy as np
import pytest

r""" Hyperlinks to mantid algorithms
LoadEmptyInstrument <https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html>
"""
from mantid.kernel import V3D

r"""
Hyperlinks to drtsans functions
ElementComponentInfo, PixelInfo, TubeInfo, TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.tubecollection import ElementComponentInfo, PixelInfo, TubeInfo, TubeCollection


@pytest.fixture(scope='function')
def simple_panel(workspace_with_instrument):
    workspace = workspace_with_instrument(axis_values=[1.0, 2.0],
                                          intensities=np.random.rand(3072).reshape((3072, 1)))
    return dict(workspace=workspace, component_info=workspace.componentInfo(), detector_info=workspace.detectorInfo())


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 3, 'n_pixels': 3,
                           'diameter': 1.0, 'height': 1.0, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
@pytest.mark.usefixtures('simple_panel')
class TestElementComponentInfo(object):

    def test_init(self, simple_panel):
        assert ElementComponentInfo(simple_panel['component_info'], 4)

    def test_decorate_component_method(self, simple_panel):
        el = ElementComponentInfo(simple_panel['component_info'], 4)
        assert el.isDetector is True  # isDetector's arity is one
        position = list(el.position)
        assert position == pytest.approx([0.0, -0.5, 0.0])
        el.setPosition(V3D(0.0, 0.0, 0.0))  # setPosition's arity is two
        assert el.position == pytest.approx([0.0, 0.0, 0.0])


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 3, 'n_pixels': 3,
                           'diameter': 1.0, 'height': 1.0, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
@pytest.mark.usefixtures('simple_panel')
class TestPixelInfo(object):

    def test_init(self, simple_panel):
        assert PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])

    def test_detector_info(self, simple_panel):
        pixel = PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])
        assert pixel.detector_info is simple_panel['detector_info']
        with pytest.raises(AttributeError):
            pixel.detector_info = 'bananas'
        # test accessibility of detectorInfo's methods
        assert pixel.l2 == pytest.approx(0.5, abs=0.01)
        assert pixel.twoTheta == pytest.approx(1.57, abs=0.01)

    def test_position(self, simple_panel):
        pixel = PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])
        assert pixel.position == pytest.approx([0, -0.5, 0])
        pixel.position = ('y', 0.0)
        assert pixel.position == pytest.approx([0, 0, 0])
        pixel.position = (0.1, 0.2, 0.3)
        assert pixel.position == pytest.approx([0.1, 0.2, 0.3])

    def test_width(self, simple_panel):
        pixel = PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])
        assert pixel.width == pytest.approx(1.0)
        pixel.width = 1.42 * pixel.width
        assert pixel.width == pytest.approx(1.42)

    def test_height(self, simple_panel):
        pixel = PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])
        assert pixel.height == pytest.approx(1.0)
        pixel.height = 1.42 * pixel.height
        assert pixel.height == pytest.approx(1.42)

    def test_area(self, simple_panel):
        pixel = PixelInfo(simple_panel['component_info'], 4, simple_panel['detector_info'])
        pixel.height, pixel.width = 1.42, 1.42
        assert pixel.area == 1.42**2


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 3, 'n_pixels': 3,
                           'diameter': 1.0, 'height': 1.0, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
@pytest.mark.usefixtures('simple_panel')
class TestTubeInfo(object):

    def test_is_valid_tube(self, simple_panel):
        assert TubeInfo(simple_panel['component_info'], 11)  # this is the index for the first tube
        with pytest.raises(ValueError):
            TubeInfo(simple_panel['component_info'], 4)  # this is a detector
        with pytest.raises(ValueError):
            TubeInfo(simple_panel['component_info'], 14)  # this is the whole panel

    def test_pixels(self, simple_panel):
        tube = TubeInfo(simple_panel['component_info'], 11)
        pixels = tube.pixels
        assert pixels[0].position == pytest.approx([1, -1.5, 0])
        assert pixels[2].position == pytest.approx([1, 0.5, 0])

    def test_getitem(self, simple_panel):
        tube = TubeInfo(simple_panel['component_info'], 11)
        assert tube[0].position == pytest.approx([1, -1.5, 0])
        assert tube[2].position == pytest.approx([1, 0.5, 0])


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 3, 'n_pixels': 3,
                           'diameter': 1.0, 'height': 1.0, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
@pytest.mark.usefixtures('simple_panel')
class TestTubeCollection(object):

    def test_init(self, simple_panel):
        collection = TubeCollection(simple_panel['workspace'], 'detector1')
        assert collection

    def test_tubes(self, simple_panel):
        collection = TubeCollection(simple_panel['workspace'], 'detector1')
        assert len(collection.tubes) == 3

    def test_getitem(self, simple_panel):
        collection = TubeCollection(simple_panel['workspace'], 'detector1')
        assert collection[0][0].position == pytest.approx([1., -1.5, 0.])  # first pixel in first tube
        assert collection[2][2].position == pytest.approx([-1., 0.5, 0.])  # last pixel in last tube

    def test_sorted(self, simple_panel):
        collection = TubeCollection(simple_panel['workspace'], 'detector1')
        sorted_tubes = collection.sorted(view='decreasing X')
        assert sorted_tubes == collection.tubes  # true for the simple_panel instrument


if __name__ == '__main__':
    pytest.main([__file__])
