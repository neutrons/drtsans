from os.path import join as pjn
import pytest
from pytest import approx
from mantid.simpleapi import Load
from drtsans.settings import (amend_config, unique_workspace_dundername as uwd)
import drtsans.sns.eqsans.dark_current as dkc


@pytest.fixture(scope='module')
def wss(reference_dir):
    with amend_config(data_dir=reference_dir.new.eqsans):
        name = pjn(reference_dir.new.eqsans, 'test_dark_current', 'data.nxs')
        # data is a Workspace2D in wavelength
        data = Load(name, OutputWorkspace=uwd())
        # dark is an EventsWorkspace in time-of-flight
        dark = Load('EQSANS_89157', OutputWorkspace=uwd())
        return dict(data=data, dark=dark)


def test_duration(wss):
    for lk in ('start_time', 'proton_charge'):
        assert dkc.duration(wss['dark'], lk).value == approx(7200.0, abs=1.0)
    with pytest.raises(AttributeError):
        dkc.duration(wss['dark'], log_key='timer')


def test_counts_in_detector(wss):
    y, e = dkc.counts_in_detector(wss['dark'])
    assert max(y) == 67.0
    assert (y[0], e[0]) == (0.0, 1.0)


if __name__ == '__main__':
    pytest.main()
