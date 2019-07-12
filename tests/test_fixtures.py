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


@pytest.mark.skip(reason="only for debugging")
@pytest.mark.parametrize('generate_sans_generic_IDF', [{'Nx':4,'Ny':4}], indirect=True)
def test_generate_IDF(generate_sans_generic_IDF):
    assert generate_sans_generic_IDF == '<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n<instrument name="GenericSANS" valid-from   ="1900-01-31 23:59:59"\n                               valid-to     ="2100-12-31 23:59:59"\n                               last-modified="2019-07-12 00:00:00">\n    <!--DEFAULTS-->\n    <defaults>\n        <length unit="metre"/>\n        <angle unit="degree"/>\n        <reference-frame>\n        <along-beam axis="z"/>\n        <pointing-up axis="y"/>\n        <handedness val="right"/>\n        <theta-sign axis="x"/>\n        </reference-frame>\n    </defaults>\n\n    <!--SOURCE-->\n    <component type="moderator">\n        <location z="-11.0"/>\n    </component>\n    <type name="moderator" is="Source"/>\n\n    <!--SAMPLE-->\n    <component type="sample-position">\n        <location y="0.0" x="0.0" z="0.0"/>\n    </component>\n    <type name="sample-position" is="SamplePos"/>\n\n    <!--RectangularDetector-->\t                       \n    <component type="panel" idstart="0" idfillbyfirst="y" idstepbyrow="4">\n        <location x="0.0" y="0.0" z="5.0" name="bank1" rot="0.0" axis-x="0" axis-y="1" axis-z="0">\n        </location>\n    </component> \n\n    <!-- Rectangular Detector Panel -->\n    <type name="panel" is="rectangular_detector" type="pixel" \n        xpixels="4" xstart="-1.5" xstep="+1.0"\n        ypixels="4" ystart="-1.5" ystep="+1.0" >\n        <properties/>\n    </type>\n\n    <!-- Pixel for Detectors-->\n    <type is="detector" name="pixel">\n        <cuboid id="pixel-shape">\n            <left-front-bottom-point y="-0.5" x="-0.5" z="0.0"/>\n            <left-front-top-point y="0.5" x="-0.5" z="0.0"/>\n            <left-back-bottom-point y="-0.5" x="-0.5" z="-0.0001"/>\n            <right-front-bottom-point y="-0.5" x="0.5" z="0.0"/>\n        </cuboid>\n        <algebra val="pixel-shape"/>\n    </type>\n\n</instrument>'


if __name__ == '__main__':
    pytest.main()
