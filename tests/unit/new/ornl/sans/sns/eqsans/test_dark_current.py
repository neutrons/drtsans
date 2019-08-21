from __future__ import (absolute_import, division, print_function)

from os.path import join as pjn
import pytest
from pytest import approx
from mantid.simpleapi import Load, SumSpectra, LoadNexus, CompareWorkspaces

from ornl.settings import (amend_config, unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
import ornl.sans.sns.eqsans.dark_current as dkc


@pytest.fixture(scope='module')
def wss(refd):
    with amend_config(data_dir=refd.new.eqsans):
        name = pjn(refd.new.eqsans, 'test_dark_current', 'data.nxs')
        # data is a Workspace2D in wavelength
        data = Load(name, OutputWorkspace=uwd())
        # dark is an EventsWorkspace in time-of-flight
        dark = Load('EQSANS_89157', OutputWorkspace=uwd())
        return dict(data=data, dark=dark)


def test_normalise_to_workspace(wss, refd):
    _w0 = dkc.normalise_to_workspace(wss['dark'], wss['data'],
                                     output_workspace=uwd())
    _w1 = SumSpectra(_w0, OutputWorkspace=uwd())
    name = pjn(refd.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=uwd())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2)]


def test_subtract_normalised_dark(wss, refd):
    name = pjn(refd.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _dark_normalised = LoadNexus(name, OutputWorkspace=uwd())
    _w0 = dkc.subtract_normalised_dark_current(wss['data'], _dark_normalised,
                                               output_workspace=uwd())
    assert SampleLogs(_w0).normalizing_duration.value == 'duration'
    _w1 = SumSpectra(_w0, OutputWorkspace=uwd())
    name = pjn(refd.new.eqsans, 'test_dark_current', 'data_minus_dark.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=uwd())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2, _dark_normalised)]


if __name__ == '__main__':
    pytest.main()
