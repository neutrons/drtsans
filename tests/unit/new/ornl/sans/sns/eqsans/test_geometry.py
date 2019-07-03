import numpy as np
import pytest
from pytest import approx
from mantid.simpleapi import Load
from ornl.sans.sns.eqsans.geometry import detector_name, center_detector


def test_center_detector(eqsans_f):
    w = Load(eqsans_f['data'])
    i = w.getInstrument()
    r_0 = np.array(i.getComponentByName(detector_name).getPos())
    r = center_detector(w, x=16, y=21, units='mm', relative=True)
    assert r_0 + np.array([.016, 0.021, 0]) == approx(r)
    r = center_detector(w, x=16, y=21, units='mm')
    assert np.array([.016, 0.021, 0]) == approx(r)


if __name__ == '__main__':
    pytest.main()
