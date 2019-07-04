import numpy as np
import pytest
from pytest import approx
from mantid.simpleapi import (LoadEventNexus, MoveInstrumentComponent)
from ornl.sans.geometry import detector_name
from ornl.sans.sns.eqsans.geometry import (center_detector,
                                           sample_aperture_diameter,
                                           source_aperture_diameter)
from ornl.sans.samplelogs import SampleLogs


def test_center_detector(eqsans_f):
    w = LoadEventNexus(eqsans_f['data'])
    r_0 = [0.1, 0.2, 0.3]
    MoveInstrumentComponent(w, X=r_0[0], Y=r_0[1], Z=r_0[2],
                            ComponentName=detector_name(w),
                            RelativePosition=False)
    i = w.getInstrument()
    r = np.array(i.getComponentByName(detector_name(i)).getPos())
    assert r == approx(r_0)
    r = center_detector(w, x=16, y=21, units='mm', relative=True)
    assert r == approx(r_0 + np.array([.016, 0.021, 0]))
    r = center_detector(w, x=16, y=21, units='mm')
    assert r == approx(np.array([.016, 0.021, r_0[2]]))
    sdd = SampleLogs(w)['sample-detector-distance'].value
    assert sdd == approx(1e3 * np.linalg.norm(r))


def test_sample_aperture_diameter():
    ws = LoadEventNexus('EQSANS_92353')
    sad = sample_aperture_diameter(ws)
    assert sad == approx(10)
    sad = SampleLogs(ws).single_value('sample-aperture-diameter')
    assert sad == approx(10)


def test_source_aperture_diameter():
    ws = LoadEventNexus('EQSANS_92353')
    sad = source_aperture_diameter(ws)
    assert sad == approx(20)
    sad = SampleLogs(ws).single_value('source-aperture-diameter')
    assert sad == approx(20)


if __name__ == '__main__':
    pytest.main()
