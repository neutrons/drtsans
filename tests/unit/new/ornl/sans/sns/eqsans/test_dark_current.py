from __future__ import (absolute_import, division, print_function)

from os.path import join as pjn
import pytest
from pytest import approx
from mantid.simpleapi import Load

import ornl.sans.sns.eqsans.dark_current as dkc


@pytest.fixture(scope='module')
def wss(refd):
    name = pjn(refd.new.eqsans, 'test_dark_current', 'EQSANS_92353.nxs')
    data = Load(name, OutputWorkspace='data')  # Matrix2DWorkspace in waveleng
    dark = Load('EQSANS_89157')  # Events workspace
    return dict(data=data, dark=dark)


def test_compute_log_ratio(wss):
    s, d = wss['data'], wss['dark']
    for lk in ('duration', 'proton_charge'):
        assert dkc.compute_log_ratio(s, d, lk) == approx(0.014, abs=0.001)
    with pytest.raises(AttributeError):
        dkc.compute_log_ratio(s, d, 'non_existent_log')


def test_duration_ratio(wss):
    s, d = wss['data'], wss['dark']
    for lk in ('duration', 'proton_charge'):
        assert dkc.duration_ratio(s, d, lk) == approx(0.014, abs=0.001)
    assert dkc.duration_ratio(s, d, 'non_existent_log') == approx(1.0)


def test_normalise_to_workspace(wss):
    s, d = wss['data'], wss['dark']
    dkc.normalise_to_workspace(d, s, 'darkn')





if __name__ == '__main__':
    pytest.main()
