from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.kernel import V3D
import pytest
from pytest import approx


@pytest.mark.skip(reason="only for debugging")
def test_eqsans_w(eqsans_w):
    pass


@pytest.mark.skip(reason="only for debugging")
def test_porasil_slice1m(porasil_slice1m):
    w = porasil_slice1m.w
    for k in w.keys():
        assert w._w[k].name() == '_porasil_slice1m_' + k
        assert w[k].name() == 'porasil_slice1m_' + k


@pytest.mark.skip(reason="only for debugging")
def test_frame_skipper(frame_skipper):
    w = frame_skipper.w
    for k in w.keys():
        assert w._w[k].name() == '_frame_skipper_' + k
        assert w[k].name() == 'frame_skipper_' + k


@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4}],
                         indirect=True)
def test_generate_IDF(generic_IDF):
    expected = '''<?xml version=\'1.0\' encoding=\'UTF-8\'?>
<instrument name="GenericSANS" valid-from   ="1900-01-31 23:59:59"
                               valid-to     ="2100-12-31 23:59:59"
                               last-modified="2019-07-12 00:00:00">
    <!--DEFAULTS-->
    <defaults>
        <length unit="metre"/>
        <angle unit="degree"/>
        <reference-frame>
        <along-beam axis="z"/>
        <pointing-up axis="y"/>
        <handedness val="right"/>
        <theta-sign axis="x"/>
        </reference-frame>
    </defaults>

    <!--SOURCE-->
    <component type="moderator">
        <location z="-11.0"/>
    </component>
    <type name="moderator" is="Source"/>

    <!--SAMPLE-->
    <component type="sample-position">
        <location y="0.0" x="0.0" z="0.0"/>
    </component>
    <type name="sample-position" is="SamplePos"/>

    <!--RectangularDetector-->
    <component type="panel" idstart="0" idfillbyfirst="y" idstepbyrow="4">
        <location x="0.0" y="0.0" z="5.0"
            name="detector1"
            rot="0.0" axis-x="0" axis-y="1" axis-z="0">
        </location>
    </component>

    <!-- Rectangular Detector Panel -->
    <type name="panel" is="rectangular_detector" type="pixel"
        xpixels="4" xstart="-1.5" xstep="+1.0"
        ypixels="4" ystart="-1.5" ystep="+1.0" >
        <properties/>
    </type>

    <!-- Pixel for Detectors-->
    <type is="detector" name="pixel">
        <cuboid id="pixel-shape">
            <left-front-bottom-point y="-0.5" x="-0.5" z="0.0"/>
            <left-front-top-point y="0.5" x="-0.5" z="0.0"/>
            <left-back-bottom-point y="-0.5" x="-0.5" z="-0.0001"/>
            <right-front-bottom-point y="-0.5" x="0.5" z="0.0"/>
        </cuboid>
        <algebra val="pixel-shape"/>
    </type>

    <parameter name="x-pixel-size">
        <value val="1000.0"/>
    </parameter>

    <parameter name="y-pixel-size">
        <value val="1000.0"/>
    </parameter>
</instrument>'''

    assert generic_IDF == expected


@pytest.mark.parametrize('generic_IDF',
                         [{}],
                         indirect=True)
def test_generate_IDF_defaults(generic_IDF):
    expected = '''<?xml version=\'1.0\' encoding=\'UTF-8\'?>
<instrument name="GenericSANS" valid-from   ="1900-01-31 23:59:59"
                               valid-to     ="2100-12-31 23:59:59"
                               last-modified="2019-07-12 00:00:00">
    <!--DEFAULTS-->
    <defaults>
        <length unit="metre"/>
        <angle unit="degree"/>
        <reference-frame>
        <along-beam axis="z"/>
        <pointing-up axis="y"/>
        <handedness val="right"/>
        <theta-sign axis="x"/>
        </reference-frame>
    </defaults>

    <!--SOURCE-->
    <component type="moderator">
        <location z="-11.0"/>
    </component>
    <type name="moderator" is="Source"/>

    <!--SAMPLE-->
    <component type="sample-position">
        <location y="0.0" x="0.0" z="0.0"/>
    </component>
    <type name="sample-position" is="SamplePos"/>

    <!--RectangularDetector-->
    <component type="panel" idstart="0" idfillbyfirst="y" idstepbyrow="3">
        <location x="0.0" y="0.0" z="5.0"
            name="detector1"
            rot="0.0" axis-x="0" axis-y="1" axis-z="0">
        </location>
    </component>

    <!-- Rectangular Detector Panel -->
    <type name="panel" is="rectangular_detector" type="pixel"
        xpixels="3" xstart="-1.0" xstep="+1.0"
        ypixels="3" ystart="-1.0" ystep="+1.0" >
        <properties/>
    </type>

    <!-- Pixel for Detectors-->
    <type is="detector" name="pixel">
        <cuboid id="pixel-shape">
            <left-front-bottom-point y="-0.5" x="-0.5" z="0.0"/>
            <left-front-top-point y="0.5" x="-0.5" z="0.0"/>
            <left-back-bottom-point y="-0.5" x="-0.5" z="-0.0001"/>
            <right-front-bottom-point y="-0.5" x="0.5" z="0.0"/>
        </cuboid>
        <algebra val="pixel-shape"/>
    </type>

    <parameter name="x-pixel-size">
        <value val="1000.0"/>
    </parameter>

    <parameter name="y-pixel-size">
        <value val="1000.0"/>
    </parameter>
</instrument>'''

    assert generic_IDF == expected


def test_generate_IDF_minimal(generic_IDF):
    assert generic_IDF


def test_generate_workspace_defaults(generic_workspace):
    ws = generic_workspace  # give it a friendly name
    assert ws
    assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
    assert ws.getNumberHistograms() == 9
    for i in range(ws.getNumberHistograms()):
        assert ws.readX(i).tolist() == [0.]
        assert ws.readY(i).tolist() == [0.]
        assert ws.readE(i).tolist() == [1.]  # SANS default


@pytest.mark.parametrize('generic_workspace',
                         [{'Nx': 2, 'Ny': 3}],
                         indirect=True)
def test_generate_workspace_mono_no_data(generic_workspace):
    ws = generic_workspace  # give it a friendly name
    assert ws
    assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
    assert ws.getNumberHistograms() == 6
    for i in range(ws.getNumberHistograms()):
        assert ws.readX(i).tolist() == [0.]
        assert ws.readY(i).tolist() == [0.]
        assert ws.readE(i).tolist() == [1.]  # SANS default


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': [42.],
                           'intensities': [[1., 4.], [9., 16.], [25., 36.]]}],
                         indirect=True)
def test_generate_workspace_monochromatic(generic_workspace):
    ws = generic_workspace  # give it a friendly name
    assert ws
    assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
    assert ws.getNumberHistograms() == 6
    for i in range(ws.getNumberHistograms()):
        assert ws.readX(i).tolist() == [42.]
    # supplied y-values
    assert ws.extractY().ravel().tolist() == [1., 4., 9., 16., 25., 36.]
    # e-values is sqrt of y
    assert ws.extractE().ravel().tolist() == [1., 2., 3., 4., 5., 6.]

    # verify particular pixels
    assert ws.readY(1) == 4.
    assert ws.readY(3) == 16.
    specInfo = ws.spectrumInfo()
    assert specInfo.position(0) == V3D(-1., -.5, 5.)  # row=0, col=0
    assert specInfo.position(3) == V3D(0., .5, 5.)    # row=1, col=0


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_units': 'tof',
                           'axis_values': [100., 8000., 16000.],
                           'intensities': [[[1., 1.], [4., 4.]], [[9., 9.], [16., 16.]], [[25., 25.], [36., 36.]]]}],
                         indirect=True)
def test_generate_workspace_tof(generic_workspace):
    ws = generic_workspace  # give it a friendly name
    assert ws
    assert ws.getAxis(0).getUnit().caption() == 'Time-of-flight'
    assert ws.getNumberHistograms() == 6
    for i in range(ws.getNumberHistograms()):
        assert ws.readX(i).tolist() == [100., 8000., 16000.]
    # supplied y-values
    assert ws.extractY().ravel().tolist() == [1., 1., 4., 4., 9., 9., 16., 16., 25., 25., 36., 36.]
    # e-values is sqrt of y
    assert ws.extractE().ravel().tolist() == [1., 1., 2., 2., 3., 3., 4., 4., 5., 5., 6., 6.]

    # verify particular pixels
    assert ws.readY(1).tolist() == [4., 4.]
    assert ws.readY(3).tolist() == [16., 16.]
    specInfo = ws.spectrumInfo()
    assert specInfo.position(0) == V3D(-1., -.5, 5.)  # row=0, col=0
    assert specInfo.position(3) == V3D(0., .5, 5.)    # row=1, col=0


def test_serve_events_workspace(serve_events_workspace):
    w1 = serve_events_workspace('EQSANS_92353')
    w2 = serve_events_workspace('EQSANS_92353')
    assert w1.name() != w2.name()
    originals = [w.name() for w in serve_events_workspace._cache.values()]
    assert w1.name() not in originals
    assert w2.name() not in originals


class TestWorkspaceWithInstrument(object):

    def test_workspace_with_instrument_defaults(self, workspace_with_instrument):
        ws = workspace_with_instrument()
        assert ws
        assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
        assert ws.getNumberHistograms() == 9
        for i in range(ws.getNumberHistograms()):
            assert ws.readX(i).tolist() == [0.]
            assert ws.readY(i).tolist() == [0.]
            assert ws.readE(i).tolist() == [1.]  # SANS default

        x, y = np.arange(9).reshape((3, 3)), np.abs(np.random.random((3, 3)))
        ws2 = workspace_with_instrument(axis_values=x, intensities=y)
        assert ws != ws2
        x, y = x.flatten(), y.flatten()
        for i in range(ws.getNumberHistograms()):
            assert ws2.readX(i).tolist() == x[i]
            assert ws2.readY(i).tolist() == y[i]
            assert ws2.readE(i).tolist() == np.sqrt(y[i])

    @pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 3, 'Ny': 2}], indirect=True)
    def test_generate_workspace_mono(self, workspace_with_instrument):
        ws = workspace_with_instrument()  # give it a friendly name
        assert ws
        assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
        assert ws.getNumberHistograms() == 6
        for i in range(ws.getNumberHistograms()):
            assert ws.readX(i).tolist() == [0.]
            assert ws.readY(i).tolist() == [0.]
            assert ws.readE(i).tolist() == [1.]  # SANS default

        ws2 = workspace_with_instrument(axis_values=[42.], view='pixel', intensities=[[1., 4.], [9., 16.], [25., 36.]])
        assert ws2
        assert ws2.getAxis(0).getUnit().caption() == 'Wavelength'
        assert ws2.getNumberHistograms() == 6
        for i in range(ws2.getNumberHistograms()):
            assert ws2.readX(i).tolist() == [42.]
        # supplied y-values
        assert ws2.extractY().ravel() == approx([1., 4., 9., 16., 25., 36.])
        # e-values is sqrt of y
        assert ws2.extractE().ravel() == approx([1., 2., 3., 4., 5., 6.])

        # verify particular pixels
        assert ws2.readY(1) == 4.
        assert ws2.readY(3) == 16.
        spectum_info = ws.spectrumInfo()
        assert spectum_info.position(0) == V3D(-1., -.5, 5.)  # row=0, col=0
        assert spectum_info.position(3) == V3D(0., .5, 5.)  # row=1, col=0

    @pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 3, 'Ny': 2}], indirect=True)
    def test_generate_workspace_tof(self, workspace_with_instrument):
        ws = workspace_with_instrument(axis_units='tof', axis_values=[100., 8000., 16000.], view='pixel',
                                       intensities=[[[1., 1.], [4., 4.]], [[9., 9.], [16., 16.]],
                                                    [[25., 25.], [36., 36.]]])
        assert ws
        assert ws.getAxis(0).getUnit().caption() == 'Time-of-flight'
        assert ws.getNumberHistograms() == 6
        for i in range(ws.getNumberHistograms()):
            assert ws.readX(i).tolist() == [100., 8000., 16000.]
        # supplied y-values
        assert ws.extractY().ravel() == approx([1., 1., 4., 4., 9., 9., 16., 16., 25., 25., 36., 36.])
        # e-values is sqrt of y
        assert ws.extractE().ravel() == approx([1., 1., 2., 2., 3., 3., 4., 4., 5., 5., 6., 6.])

        # verify particular pixels
        assert ws.readY(1).tolist() == [4., 4.]
        assert ws.readY(3).tolist() == [16., 16.]
        spectrum_info = ws.spectrumInfo()
        assert spectrum_info.position(0) == V3D(-1., -.5, 5.)  # row=0, col=0
        assert spectrum_info.position(3) == V3D(0., .5, 5.)  # row=1, col=0

    @pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 3, 'Ny': 4}], indirect=True)
    def test_correct_pixel_id(self, workspace_with_instrument):
        intensities = np.array([[3,  7, 11],
                                [2,  6, 10],
                                [1,  5,  9],
                                [0,  4,  8]])  # three tubes, intensities are pixel id's
        ws = workspace_with_instrument(axis_values=[1.], intensities=intensities, view='array')
        assert ws.extractY().ravel().tolist() == list(range(12))
        intensities = np.array([[0, 1, 2, 3],
                                [4, 5, 6, 7],
                                [8, 9, 10, 11]])  # pixel view, each tube along the column axis
        ws = workspace_with_instrument(axis_values=[1.], intensities=intensities, view='pixel')
        assert ws.extractY().ravel().tolist() == list(range(12))

    @pytest.mark.parametrize('workspace_with_instrument', [dict(name='Theod-osius_IV4')], indirect=True)
    def test_instrument_name(self, workspace_with_instrument):
        ws = workspace_with_instrument()
        assert ws.getInstrument().getName() == 'Theod-osius_IV4'


if __name__ == '__main__':
    pytest.main()
