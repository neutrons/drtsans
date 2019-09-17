import pytest
from pytest import approx
import numpy as np

from mantid.simpleapi import (Rebin, SumSpectra)
from drtsans.settings import (amend_config, unique_workspace_name as uwn)
from drtsans.tof.eqsans.load import load_events, load_events_monitor


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


if __name__ == '__main__':
    pytest.main()
