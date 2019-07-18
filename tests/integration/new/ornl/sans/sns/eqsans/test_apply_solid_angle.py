import os
import pytest
from ornl.sans.sns import eqsans
from mantid.simpleapi import LoadEmptyInstrument, SolidAngle
from mantid.kernel import V3D
import numpy as np
from copy import deepcopy


@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 3, 'Ny': 3, 'dx': 0.00425,
                           'dy': 0.0055, 'xr': 0.32, 'yr': -0.16}],
                         indirect=True)
def test_solid_angle(generate_sans_generic_IDF):
    tmp = open(r'/tmp/GenericSANS_Definition.xml', 'w')
    tmp.write(generate_sans_generic_IDF)
    tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS')
    ws.dataY(4)[0] = 156.
    ws.dataE(4)[0] = np.sqrt(156.)

    assert ws.dataY(4)[0] == 156.
    assert ws.dataE(4)[0] == np.sqrt(156.)

    d_info = ws.detectorInfo()
    s_info = ws.spectrumInfo()

    # c = [0., 0., 5.000]
    s = s_info.samplePosition()
    assert s == V3D(0., 0., 0.)
    r = d_info.position(4)
    assert r == V3D(0.320, -0.160, 5.)
    print(r.norm())
    cos_two_theta = np.cos(d_info.twoTheta(4))
    assert cos_two_theta == pytest.approx(0.9974497886)

    b = deepcopy(r)
    b[1] = 0.0
    cos_alpha = b.scalar_prod(r)/(b.norm()*r.norm())
    assert cos_alpha == pytest.approx(0.9994904782)

    ws2 = SolidAngle(InputWorkspace=str(ws), Method='VerticalTube')
    assert ws2.dataY(4)[0] == pytest.approx(9.29763209E-07)

    ws = eqsans.apply_solid_angle_correction(str(ws))

    assert ws.dataY(4)[0] == pytest.approx(167784655.715)
    assert ws.dataE(4)[0] == pytest.approx(13433523.5782)

    os.remove(r'/tmp/GenericSANS_Definition.xml')


if __name__ == '__main__':
    pytest.main()
