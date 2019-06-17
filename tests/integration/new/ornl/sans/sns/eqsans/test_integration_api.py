from __future__ import (absolute_import, division, print_function)

import pytest
from mantid.dataobjects import EventWorkspace
from ornl.sans.sns import eqsans
from ornl.settings import unique_workspace_name as uwn


#               run        num-events max-tof #SDD
run_sets = (('EQSANS_92353', 262291, 33332.5, 4000.0),)  # old IDF,
"""To-Do after PR new IDF accepted in mantid master
run_sets = (('EQSANS_92353', 262291, 33332.5, 0.0),  # old IDF, no translation
            ('EQSANS_102616', 0.0, 0.0, 0.0))  # new IDF, translation
"""


@pytest.mark.parametrize('run_set', run_sets)
def test_load_events(run_set):
    ws = eqsans.load_events(run_set[0], output_workspace=uwn())
    assert isinstance(ws, EventWorkspace)
    assert ws.getNumberEvents() == run_set[1]
    assert ws.getTofMax() == run_set[2]
    # assert distance of detector1 same as that in detectorZ of the logs
    instrument = ws.getInstrument()
    det = instrument.getComponentByName(eqsans.detector_name)
    d1 = det.getDistance(instrument.getSample())
    assert run_set[3] == pytest.approx(d1 * 1000)


if __name__ == '__main__':
    pytest.main()
