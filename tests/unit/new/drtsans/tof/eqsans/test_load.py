import pytest
from pytest import approx
import numpy as np

from mantid.simpleapi import (Rebin, SumSpectra, mtd)
from drtsans.settings import (amend_config, unique_workspace_name as uwn)
from drtsans.tof.eqsans.load import load_events, load_events_monitor, merge_data
from drtsans.tof.eqsans.correct_frame import transform_to_wavelength
from drtsans.samplelogs import SampleLogs


def test_load_events(reference_dir):
    # default workspace name is file hint
    with amend_config(data_dir=reference_dir.new.eqsans):
        ws_test_load_events = load_events('EQSANS_92353')
    assert ws_test_load_events.name() == 'EQSANS_92353'

    ws_name = uwn()
    with amend_config(data_dir=reference_dir.new.eqsans):
        ws = load_events('EQSANS_92353', output_workspace=ws_name)
    assert ws.name() == ws_name

    assert ws.getTofMin() == pytest.approx(11410, abs=1)
    assert ws.getTofMax() == pytest.approx(61439, abs=1)

    ws = Rebin(ws, Params=[10000, 1000, 62000], PreserveEvents=False)
    ws = SumSpectra(ws)
    assert len(np.nonzero(ws.dataY(0))[0]) == 36


def test_load_events_monitor(reference_dir):
    # Raises for a run in skip frame mode
    with pytest.raises(RuntimeError, match='cannot correct monitor'):
        load_events_monitor('EQSANS_92353', data_dir=reference_dir.new.eqsans)

    w = load_events_monitor('EQSANS_88901', data_dir=reference_dir.new.eqsans)
    assert w.name() == 'EQSANS_88901_monitors'
    assert w.getSpectrum(0).getTofMin() == approx(30680, abs=1)
    assert w.getSpectrum(0).getTofMax() == approx(47346, abs=1)


def test_merge_Data(reference_dir):
    ws0 = load_events('EQSANS_101595', data_dir=reference_dir.new.eqsans)
    ws0 = transform_to_wavelength(ws0)
    ws1 = load_events('EQSANS_104088', data_dir=reference_dir.new.eqsans)
    ws1 = transform_to_wavelength(ws1)
    ws2 = load_events('EQSANS_105428', data_dir=reference_dir.new.eqsans)
    ws2 = transform_to_wavelength(ws2)

    sample_logs0 = SampleLogs(ws0)
    sample_logs1 = SampleLogs(ws1)
    sample_logs2 = SampleLogs(ws2)

    merged_workspaces = merge_data([ws0, ws1, ws2],  output_workspace="merged")

    merged_sample_logs = SampleLogs(merged_workspaces)

    # Check duration increase as the sum
    assert sample_logs0.duration.value == pytest.approx(215.531066895, abs=1e-9)
    assert sample_logs1.duration.value == pytest.approx(289.029266357, abs=1e-9)
    assert sample_logs2.duration.value == pytest.approx(72.3323135376, abs=1e-9)
    assert merged_sample_logs.duration.value == pytest.approx(215.531066895 + 289.029266357 + 72.3323135376, abs=1e-9)

    # Check proton charge is correct
    assert sample_logs0.getProtonCharge() == pytest.approx(83.37074628055555, abs=1e-9)
    assert sample_logs1.getProtonCharge() == pytest.approx(111.1237739861111, abs=1e-9)
    assert sample_logs2.getProtonCharge() == pytest.approx(27.799524525, abs=1e-9)
    assert merged_sample_logs.getProtonCharge() == pytest.approx(83.37074628055555, + 111.1237739861111 + 27.799524525,
                                                                 abs=1e-9)

    # Check Time Series properties increase length
    assert sample_logs0.proton_charge.size() == 12933
    assert sample_logs1.proton_charge.size() == 17343
    assert sample_logs2.proton_charge.size() == 4341
    assert merged_sample_logs.proton_charge.size() == 12933 + 17343 + 4341

    # Check integrated intensity increases as the total sum
    assert mtd[str(ws0)].extractY().sum() == 289051
    assert mtd[str(ws1)].extractY().sum() == 1318463
    assert mtd[str(ws2)].extractY().sum() == 65694
    assert mtd[str(merged_workspaces)].extractY().sum() == 289051 + 1318463 + 65694


if __name__ == '__main__':
    pytest.main([__file__])
