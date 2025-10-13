from mantid.simpleapi import LoadEmptyInstrument, mtd
import pytest

from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans.meta_data import is_sample_run


def test_is_sample_run():
    workspace_name = mtd.unique_hidden_name()
    LoadEmptyInstrument(InstrumentName="EQSANS", OutputWorkspace=workspace_name)

    # no logs
    assert is_sample_run(workspace_name) is False

    # empty logs
    sample_logs = SampleLogs(workspace_name)
    sample_logs.insert_time_series("SampleName", elapsed_times=[0.0], values=[""])
    sample_logs.insert_time_series("SampleFormula", elapsed_times=[0.0], values=[""])
    assert sample_logs.single_value("SampleName") == ""
    assert sample_logs.single_value("SampleFormula") == ""
    assert is_sample_run(workspace_name) is False

    # only one non-empty log
    sample_logs.insert_time_series("SampleName", elapsed_times=[0.0], values=["name"])  # overwrite log entry
    assert sample_logs.single_value("SampleName") == "name"
    assert is_sample_run(workspace_name) is False

    # both logs are non-empty
    sample_logs.insert_time_series("SampleFormula", elapsed_times=[0.0], values=["formula"])  # overwrite log entry
    assert sample_logs.single_value("SampleFormula") == "formula"
    assert is_sample_run(workspace_name) is True


if __name__ == "__main__":
    pytest.main([__file__])
