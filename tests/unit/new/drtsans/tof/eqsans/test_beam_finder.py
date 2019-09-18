import numpy as np
import pytest
from pytest import approx
from mantid.simpleapi import (LoadEventNexus, MoveInstrumentComponent)

from drtsans.settings import amend_config, unique_workspace_name as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import detector_name
from drtsans.tof.eqsans.beam_finder import center_detector, find_beam_center


def test_find_beam_center(reference_dir):
    with amend_config(data_dir=reference_dir.new.eqsans):
        w = LoadEventNexus(Filename='EQSANS_92160', OutputWorkspace=uwd())
    assert find_beam_center(w) == approx((0.025, 0.013), abs=1e-3)


def test_center_detector(reference_dir):
    #
    # Translate detector according to method center_of_mass
    #
    with amend_config(data_dir=reference_dir.new.eqsans):
        w = LoadEventNexus(Filename='EQSANS_92160', OutputWorkspace=uwd())
    r = center_detector(w, method='center_of_mass')
    assert r == approx((-0.0255, -0.0131, 0), abs=0.0001)
    #
    # Translate detector according to method center_of_mass plus additional
    # translation
    #
    MoveInstrumentComponent(w, X=0, Y=0, Z=0, ComponentName=detector_name(w),
                            RelativePosition=False)
    r = center_detector(w, method='center_of_mass', x=0.01, y=0.02,
                        relative=True)
    assert r == approx((-0.0155, 0.0069, 0), abs=0.001)
    r_0 = [0.1, 0.2, 0.3]
    MoveInstrumentComponent(w, X=r_0[0], Y=r_0[1], Z=r_0[2],
                            ComponentName=detector_name(w),
                            RelativePosition=False)
    #
    # Translate detector by XY amount
    #
    r = center_detector(w, method=None, x=16, y=21, unit='mm', relative=True)
    assert r == approx(r_0 + np.array([.016, 0.021, 0]))
    #
    # Translate detector to XY position, leave Z unchanged
    #
    r = center_detector(w, x=16, y=21, unit='mm')
    assert r == approx(np.array([.016, 0.021, r_0[2]]))
    sdd = SampleLogs(w)['sample-detector-distance'].value
    assert sdd == approx(1e3 * np.linalg.norm(r))


if __name__ == '__main__':
    pytest.main()