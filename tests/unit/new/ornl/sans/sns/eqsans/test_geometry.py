import numpy as np
import pytest
from pytest import approx
from mantid.simpleapi import Load
from ornl.sans.geometry import detector_name
from ornl.sans.sns.eqsans.geometry import center_detector
from ornl.sans.samplelogs import SampleLogs


def test_center_detector(eqsans_f):
    w = Load(eqsans_f['data'])
    i = w.getInstrument()
    r_0 = np.array(i.getComponentByName(detector_name(i)).getPos())
    r = center_detector(w, x=16, y=21, units='mm', relative=True)
    assert r_0 + np.array([.016, 0.021, 0]) == approx(r)
    r = center_detector(w, x=16, y=21, units='mm')
    assert np.array([.016, 0.021, 0]) == approx(r)
    sdd = SampleLogs(w)['sample-detector-distance'].value
    assert sdd == approx(1e3 * np.linalg.norm(np.array([.016, 0.021, 0])))


if __name__ == '__main__':
    pytest.main()
