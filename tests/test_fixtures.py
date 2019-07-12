from __future__ import (absolute_import, division, print_function)

import pytest


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


@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 4, 'Ny': 4}],
                         indirect=True)
def test_generate_IDF(generate_sans_generic_IDF):
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
            name="bank1"
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

</instrument>'''

    assert generate_sans_generic_IDF == expected


if __name__ == '__main__':
    pytest.main()
