from mantid.kernel import DateAndTime, FloatTimeSeriesProperty
from mantid.simpleapi import CreateWorkspace
import pytest

from drtsans.tof.eqsans.samplelogs import (
    DEFAULT_FREQUENCY,
    FALLBACK_FREQUENCY_LOG,
    FREQUENCY_LOG,
    get_frequency,
)


def add_frequency_log(workspace, name, values):
    log = FloatTimeSeriesProperty(name)
    for index, value in enumerate(values):
        log.addValue(DateAndTime(f"2020-01-01T00:00:0{index}"), value)
    workspace.mutableRun().addProperty(name, log, True)


def test_get_frequency_returns_primary_mean(clean_workspace):
    workspace = CreateWorkspace(DataX=[0.0, 1.0], DataY=[1.0])
    clean_workspace(workspace)
    add_frequency_log(workspace, FREQUENCY_LOG, [50.0, 70.0])
    add_frequency_log(workspace, FALLBACK_FREQUENCY_LOG, [30.0, 30.0])

    assert get_frequency(workspace) == pytest.approx(60.0)


def test_get_frequency_uses_fallback_when_primary_is_missing(clean_workspace):
    workspace = CreateWorkspace(DataX=[0.0, 1.0], DataY=[1.0])
    clean_workspace(workspace)
    add_frequency_log(workspace, FALLBACK_FREQUENCY_LOG, [50.0, 70.0])

    assert get_frequency(workspace) == pytest.approx(60.0)


def test_get_frequency_uses_fallback_when_primary_is_zero(clean_workspace):
    workspace = CreateWorkspace(DataX=[0.0, 1.0], DataY=[1.0])
    clean_workspace(workspace)
    add_frequency_log(workspace, FREQUENCY_LOG, [0.0])
    add_frequency_log(workspace, FALLBACK_FREQUENCY_LOG, [59.0, 61.0])

    assert get_frequency(workspace) == pytest.approx(60.0)


def test_get_frequency_uses_default_when_both_logs_are_unusable(clean_workspace):
    workspace = CreateWorkspace(DataX=[0.0, 1.0], DataY=[1.0])
    clean_workspace(workspace)
    add_frequency_log(workspace, FREQUENCY_LOG, [0.0])
    add_frequency_log(workspace, FALLBACK_FREQUENCY_LOG, [0.0])

    assert get_frequency(workspace) == pytest.approx(DEFAULT_FREQUENCY)
