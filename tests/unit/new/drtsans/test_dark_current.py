from os.path import join as pjn
import pytest
from pytest import approx

r""" Links to mantid algorithms
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
"""
from mantid.simpleapi import Load

r"""
Hyperlinks to drtsans functions
amend_config, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
duration, counts_in_detector <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/dark_current.py>
"""  # noqa: E501
from drtsans.settings import amend_config, unique_workspace_dundername
from drtsans.dark_current import duration, counts_in_detector


@pytest.fixture(scope='module')
def wss(reference_dir):
    with amend_config(data_dir=reference_dir.new.eqsans):
        name = pjn(reference_dir.new.eqsans, 'test_dark_current', 'data.nxs')
        # data is a Workspace2D in wavelength
        data = Load(name, OutputWorkspace=unique_workspace_dundername())
        # dark is an EventsWorkspace in time-of-flight
        dark = Load('EQSANS_89157', OutputWorkspace=unique_workspace_dundername())
        return dict(data=data, dark=dark)


def test_duration(wss):
    for lk in ('start_time', 'duration'):
        assert duration(wss['dark'], lk).value == approx(7200.0, abs=1.0)
    with pytest.raises(AttributeError):
        duration(wss['dark'], log_key='timer')


def test_counts_in_detector(wss):
    y, e = counts_in_detector(wss['dark'])
    assert max(y) == 67.0
    assert (y[0], e[0]) == (0.0, 1.0)


if __name__ == '__main__':
    pytest.main([__file__])
