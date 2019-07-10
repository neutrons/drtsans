import numpy as np
import pytest
from pytest import approx
from mantid.simpleapi import (LoadEventNexus, MoveInstrumentComponent)
from ornl.settings import (amend_config, unique_workspace_dundername as uwd)
from ornl.sans.geometry import detector_name
from ornl.sans.sns.eqsans.geometry import (center_detector, direct_beam_center,
                                           sample_aperture_diameter,
                                           source_aperture_diameter)
from ornl.sans.samplelogs import SampleLogs


def test_center_detector(eqsans_f):
    #
    # Translate detector according to method center_of_mass
    #
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        w = LoadEventNexus(Filename='EQSANS_92160', OutputWorkspace=uwd())
    r = center_detector(w, method='center_of_mass')
    assert r == approx((0.025, 0.013, 0), abs=0.001)
    #
    # Translate detector according to method center_of_mass plus additional
    # translation
    #
    MoveInstrumentComponent(w, X=0, Y=0, Z=0, ComponentName=detector_name(w),
                            RelativePosition=False)
    r = center_detector(w, method='center_of_mass',
                        x=0.01, y=0.02, relative=True)
    assert r == approx((0.035, 0.033, 0), abs=0.001)
    r_0 = [0.1, 0.2, 0.3]
    MoveInstrumentComponent(w, X=r_0[0], Y=r_0[1], Z=r_0[2],
                            ComponentName=detector_name(w),
                            RelativePosition=False)
    #
    # Translate detector by XY amount
    #
    r = center_detector(w, method=None, x=16, y=21, units='mm', relative=True)
    assert r == approx(r_0 + np.array([.016, 0.021, 0]))
    #
    # Translate detector to XY position, leave Z unchanged
    #
    r = center_detector(w, x=16, y=21, units='mm')
    assert r == approx(np.array([.016, 0.021, r_0[2]]))
    sdd = SampleLogs(w)['sample-detector-distance'].value
    assert sdd == approx(1e3 * np.linalg.norm(r))


def test_direct_beam_center():
    with amend_config({'datasearch.searcharchive': 'hfir,sns'}):
        w = LoadEventNexus(Filename='EQSANS_92160', OutputWorkspace=uwd())
    assert direct_beam_center(w) == approx((0.025, 0.013), abs=1e-3)


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
