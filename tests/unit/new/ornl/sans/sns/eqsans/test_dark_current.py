from __future__ import (absolute_import, division, print_function)

from os.path import join as pjn
import pytest
from pytest import approx
import numpy as np
from mantid.simpleapi import Load, SumSpectra, LoadNexus, CompareWorkspaces, CreateWorkspace

from ornl.settings import (amend_config, unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
import ornl.sans.sns.eqsans.dark_current as dkc


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


def test_normalise_to_workspace(wss, reference_dir):
    _w0 = dkc.normalise_to_workspace(wss['dark'], wss['data'],
                                     output_workspace=uwd())
    _w1 = SumSpectra(_w0, OutputWorkspace=uwd())
    name = pjn(reference_dir.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=uwd())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2)]


def test_subtract_normalised_dark(wss, reference_dir):
    name = pjn(reference_dir.new.eqsans, 'test_dark_current', 'dark_norm_sum.nxs')
    _dark_normalised = LoadNexus(name, OutputWorkspace=uwd())
    _w0 = dkc.subtract_normalised_dark_current(wss['data'], _dark_normalised,
                                               output_workspace=uwd())
    assert SampleLogs(_w0).normalizing_duration.value == 'duration'
    _w1 = SumSpectra(_w0, OutputWorkspace=uwd())
    name = pjn(reference_dir.new.eqsans, 'test_dark_current', 'data_minus_dark.nxs')
    _w2 = LoadNexus(name, OutputWorkspace=uwd())
    assert CompareWorkspaces(_w1, _w2)
    [_w.delete() for _w in (_w0, _w1, _w2, _dark_normalised)]


def test_flatten_TOF():
    tof = [1., 2., 3., 4.] * 9  # wavelength boundaries
    cts = [23, 5, 15, 18, 50, 13, 9, 7, 15,
           48, 41, 34, 79, 45, 33, 85, 78, 1,
           50, 20, 105, 53, 23, 45, 47, 30, 45]
    err = np.sqrt(cts)
    ws = CreateWorkspace(DataX=tof,
                         DataY=cts,
                         DataE=err,
                         NSpec=9)
    y, e = dkc.counts_in_detector(ws)
    expected_counts = [43, 81, 31, 123, 157, 164, 175, 121, 122]
    expected_errors = np.sqrt(expected_counts)
    assert np.allclose(y, expected_counts)
    assert np.allclose(e, expected_errors)


if __name__ == '__main__':
    pytest.main()
