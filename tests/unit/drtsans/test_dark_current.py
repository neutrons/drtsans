from os.path import join as pjn
import pytest
from pytest import approx

r""" Links to mantid algorithms
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
amend_config <https://docs.mantidproject.org/nightly/api/python/mantid/kernel/AmendConfig.html>
"""
from mantid.simpleapi import Load, mtd
from mantid.kernel import amend_config

r"""
Hyperlinks to drtsans functions
duration, counts_in_detector <https://github.com/neutrons/drtsans/blob/next/src/drtsans/dark_current.py>
"""  # noqa: E501

from drtsans.dark_current import duration, counts_in_detector


@pytest.fixture(scope="module")
def workspaces(datarepo_dir):
    with amend_config(data_dir=datarepo_dir.eqsans):
        name = pjn(datarepo_dir.eqsans, "test_dark_current", "data.nxs")
        # data is a Workspace2D in wavelength
        data = Load(name, OutputWorkspace=mtd.unique_hidden_name())
        # dark is an EventsWorkspace in time-of-flight
        dark = Load("EQSANS_89157", OutputWorkspace=mtd.unique_hidden_name())
        return dict(data=data, dark=dark)


@pytest.mark.datarepo
def test_duration(workspaces):
    for log_key in ("duration", "start_time"):
        dark_duration = duration(workspaces["dark"], log_key=log_key)
        assert dark_duration.value == approx(7200.0, abs=1.0)
    with pytest.raises(AttributeError):
        duration(workspaces["dark"], log_key="timer")


@pytest.mark.datarepo
def test_counts_in_detector(workspaces):
    y, e = counts_in_detector(workspaces["dark"])
    assert max(y) == 67.0
    assert (y[0], e[0]) == (0.0, 1.0)


if __name__ == "__main__":
    pytest.main([__file__])
