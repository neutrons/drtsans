import pytest

r""" Links to mantid algorithms
LoadHFIRSANS <https://docs.mantidproject.org/nightly/algorithms/LoadHFIRSANS-v1.html>
"""
from mantid.simpleapi import LoadHFIRSANS, DeleteWorkspaces
from mantid import mtd

r"""
Hyperlinks to drtsans functions
SampleLogs <https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py>
time <https://github.com/neutrons/drtsans/blob/next/src/drtsans/mono/normalization.py>
subtract_dark_current <https://github.com/neutrons/drtsans/blob/next/src/drtsans/mono/dark_current.py>
"""  # noqa: E501
from drtsans.samplelogs import SampleLogs
from drtsans.dark_current import duration
from drtsans.mono.dark_current import subtract_dark_current, normalize_dark_current


@pytest.mark.datarepo
def test_dark_current_workspace(gpsans_f):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # First read the data
    sample_workspace = mtd.unique_hidden_name()
    LoadHFIRSANS(Filename=gpsans_f["sample_scattering"], OutputWorkspace=sample_workspace)

    # second read and normalize the dark current
    dark_current_workspace = mtd.unique_hidden_name()
    LoadHFIRSANS(Filename=gpsans_f["dark_current"], OutputWorkspace=dark_current_workspace)
    normalized_dark_current = mtd.unique_hidden_name()
    normalize_dark_current(dark_current_workspace, output_workspace=normalized_dark_current)

    # third let's a DC subraction
    sample_subtracted = mtd.unique_hidden_name()
    subtract_dark_current(sample_workspace, dark_current_workspace, output_workspace=sample_subtracted)

    # Let's test:
    duration_log_key = SampleLogs(normalized_dark_current).normalizing_duration.value
    sample_duration = duration(sample_workspace, log_key=duration_log_key).value

    sample_sample_value = mtd[sample_workspace].dataY(612)[0]
    normalized_dark_current_sample_value = mtd[normalized_dark_current].dataY(612)[0]
    sample_subtracted_sample_value = mtd[sample_subtracted].dataY(612)[0]

    test_value = sample_sample_value - sample_duration * normalized_dark_current_sample_value
    assert sample_subtracted_sample_value == pytest.approx(test_value, abs=1.0e-6)
    DeleteWorkspaces(
        [
            sample_workspace,
            dark_current_workspace,
            normalized_dark_current,
            sample_subtracted,
        ]
    )


@pytest.mark.datarepo
def test_dark_current_filename(gpsans_f):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # First read the data
    sample_workspace = mtd.unique_hidden_name()
    LoadHFIRSANS(Filename=gpsans_f["sample_scattering"], OutputWorkspace=sample_workspace)

    # second read and normalize the dark current
    dark_current_workspace = mtd.unique_hidden_name()
    LoadHFIRSANS(Filename=gpsans_f["dark_current"], OutputWorkspace=dark_current_workspace)
    normalized_dark_current = mtd.unique_hidden_name()
    normalize_dark_current(dark_current_workspace, output_workspace=normalized_dark_current)

    # third let's a DC subraction
    sample_subtracted = mtd.unique_hidden_name()
    subtract_dark_current(sample_workspace, gpsans_f["dark_current"], output_workspace=sample_subtracted)

    # Let's test:
    duration_log_key = SampleLogs(normalized_dark_current).normalizing_duration.value
    sample_duration = duration(sample_workspace, log_key=duration_log_key).value

    sample_sample_value = mtd[sample_workspace].dataY(612)[0]
    normalized_dark_current_sample_value = mtd[normalized_dark_current].dataY(612)[0]
    sample_subtracted_sample_value = mtd[sample_subtracted].dataY(612)[0]

    test_value = sample_sample_value - sample_duration * normalized_dark_current_sample_value
    assert sample_subtracted_sample_value == pytest.approx(test_value, abs=1.0e-6)
    DeleteWorkspaces(
        [
            sample_workspace,
            dark_current_workspace,
            normalized_dark_current,
            sample_subtracted,
        ]
    )
