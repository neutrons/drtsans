from __future__ import (absolute_import, division, print_function)

from os.path import join as pjn
import pytest
from pytest import approx
from mantid.simpleapi import Load, SumSpectra, LoadNexus, CompareWorkspaces

from ornl.settings import amend_config, unique_workspace_name
from ornl.sans.samplelogs import SampleLogs
import ornl.sans.sns.eqsans.dark_current as dkc


@pytest.fixture(scope='module')
def wss(refd):
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        name = pjn(refd.new.eqsans, 'test_dark_current', 'data.nxs')
        data = Load(name, OutputWorkspace=unique_workspace_name())  # Matrix2DWorkspace in wavel
        dark = Load('EQSANS_89157', OutputWorkspace=unique_workspace_name())  # Events workspace
        return dict(data=data, dark=dark)


def test_duration(wss):
    assert dkc.duration(wss['dark']).value == approx(7200.0, abs=1.0)


def test_counts_in_detector(wss):
    y, e = dkc.counts_in_detector(wss['dark'])
    assert max(y) == 67.0
    assert (y[0], e[0]) == (0.0, 1.0)


def test_normalise_to_workspace(wss, refd):
    _w0 = dkc.normalise_to_workspace(wss['dark'], wss['data'],
                                     unique_workspace_name())
    _w1 = SumSpectra(_w0, OutputWorkspace=unique_workspace_name())
    assert SampleLogs(_w1).normalizing_duration.value == 'duration'
    name = pjn(refd.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=unique_workspace_name())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2)]


def test_subtract_normalized_dark(wss, refd):
    name = pjn(refd.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _dark_normalized = LoadNexus(name, OutputWorkspace=unique_workspace_name())
    _w0 = dkc.subtract_normalized_dark(wss['data'], _dark_normalized,
                                       unique_workspace_name())
    assert SampleLogs(_w0).normalizing_duration.value == 'duration'
    _w1 = SumSpectra(_w0, OutputWorkspace=unique_workspace_name())
    name = pjn(refd.new.eqsans, 'test_dark_current', 'data_minus_dark.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=unique_workspace_name())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2, _dark_normalized)]


if __name__ == '__main__':
    pytest.main()
