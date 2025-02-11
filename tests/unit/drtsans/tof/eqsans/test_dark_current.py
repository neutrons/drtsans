from os.path import join as pjn
import pytest
import numpy as np

r""" Links to mantid algorithms
CompareWorkspaces <https://docs.mantidproject.org/nightly/algorithms/CompareWorkspaces-v1.html>
CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
LoadNexus <https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html>
SumSpectra <https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html>
amend_config <https://docs.mantidproject.org/nightly/api/python/mantid/kernel/AmendConfig.html>
"""
from mantid.simpleapi import (
    CompareWorkspaces,
    CreateWorkspace,
    Load,
    LoadNexus,
    SumSpectra,
    mtd,
)
from mantid.kernel import amend_config

r"""
Hyperlinks to drtsans functions
SampleLogs <https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py>
dark_current <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof.eqsans/dark_current.py>
"""  # noqa: E501
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import dark_current


def test_flatten_TOF(clean_workspace):
    r"""
    Check that the counts are added together in each spectra

    Function tested: drtsans.tof.eqsans.dark_current.counts_in_detector
    Undelying Mantid algorithms:
        Integration https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html
        Transpose   https://docs.mantidproject.org/nightly/algorithms/Transpose-v1.html

    dev - Andrei Savici <saviciat@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>
    """
    # create the workspace
    tof = [1.0, 2.0, 3.0, 4.0] * 9  # wavelength boundaries
    cts = [
        23,
        5,
        15,
        18,
        50,
        13,
        9,
        7,
        15,
        48,
        41,
        34,
        79,
        45,
        33,
        85,
        78,
        1,
        50,
        20,
        105,
        53,
        23,
        45,
        47,
        30,
        45,
    ]
    err = np.sqrt(cts)
    ws = CreateWorkspace(DataX=tof, DataY=cts, DataE=err, NSpec=9)
    clean_workspace(ws)
    # run the function
    y, e = dark_current.counts_in_detector(ws)
    # check the results
    expected_counts = [43, 81, 31, 123, 157, 164, 175, 121, 122]
    expected_errors = np.sqrt(expected_counts)
    assert np.allclose(y, expected_counts)
    assert np.allclose(e, expected_errors)


@pytest.fixture(scope="module")
def wss(datarepo_dir):
    with amend_config(data_dir=datarepo_dir.eqsans):
        name = pjn(datarepo_dir.eqsans, "test_dark_current", "data.nxs")
        # data is a Workspace2D in wavelength
        data = Load(name, OutputWorkspace=mtd.unique_hidden_name())
        # dark is an EventsWorkspace in time-of-flight
        dark = Load("EQSANS_89157", OutputWorkspace=mtd.unique_hidden_name())
        return dict(data=data, dark=dark)


@pytest.mark.datarepo
def test_normalize_to_workspace(wss, datarepo_dir, temp_workspace_name, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    pytest.skip("This test fails, defect written up in EWM Defect 2841")
    _w0 = dark_current.normalize_dark_current(wss["dark"], wss["data"], output_workspace=temp_workspace_name())
    _w1 = SumSpectra(_w0, OutputWorkspace=temp_workspace_name())
    name = pjn(datarepo_dir.eqsans, "test_dark_current", "dark_norm_sum.nxs")
    _w2 = LoadNexus(name, OutputWorkspace=temp_workspace_name())
    result, messages = CompareWorkspaces(_w1, _w2)
    clean_workspace(messages)
    assert result


@pytest.mark.datarepo
def test_subtract_normalized_dark(wss, datarepo_dir, temp_workspace_name, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    file_path = pjn(datarepo_dir.eqsans, "test_dark_current", "dark_norm_sum.nxs")
    dark_normalized = LoadNexus(file_path, OutputWorkspace=temp_workspace_name())
    data_normalized = dark_current.subtract_normalized_dark_current(
        wss["data"], dark_normalized, output_workspace=temp_workspace_name()
    )
    assert SampleLogs(data_normalized).normalizing_duration.value == "duration"
    summed_normalized = SumSpectra(data_normalized, OutputWorkspace=temp_workspace_name())

    # Compare to stored data
    file_path = pjn(datarepo_dir.eqsans, "test_dark_current", "data_minus_dark.nxs")
    stored_summed_normalized = LoadNexus(file_path, OutputWorkspace=temp_workspace_name())
    result, messages = CompareWorkspaces(summed_normalized, stored_summed_normalized)
    clean_workspace(messages)
    assert result


if __name__ == "__main__":
    pytest.main([__file__])
