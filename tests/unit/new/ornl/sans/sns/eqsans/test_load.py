import pytest
import numpy as np

from mantid.simpleapi import (Rebin, SumSpectra)
from ornl.settings import unique_workspace_name as uwn
from ornl.sans.sns.eqsans.load import load_events


def test_load_events():
    # default workspace name is file hint
    ws_test_load_events = load_events('EQSANS_92353')
    assert ws_test_load_events.name() == 'EQSANS_92353'

    ws_name = uwn()
    ws = load_events('EQSANS_92353', output_workspace=ws_name)
    assert ws.name() == ws_name

    assert ws.getTofMin() == pytest.approx(11410, abs=1)
    assert ws.getTofMax() == pytest.approx(61412, abs=1)

    ws = Rebin(ws, Params=[10000, 1000, 62000], PreserveEvents=False)
    ws = SumSpectra(ws)
    assert len(np.nonzero(ws.dataY(0))[0]) == 36


if __name__ == '__main__':
    pytest.main()
