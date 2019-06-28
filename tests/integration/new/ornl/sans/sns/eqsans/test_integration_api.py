"""
    Test top-level API
"""
import pytest
from pytest import approx
from mantid.dataobjects import EventWorkspace

# public API
from ornl.sans.sns import eqsans

# protected API
from ornl.settings import (namedtuplefy, unique_workspace_name as uwn)
from ornl.sans.samplelogs import SampleLogs


keys = ('run', 'num_events', 'nominal_sdd', 'sdd', 'ssd', 'min_tof', 'max_tof',
        'skip_frame')
values = (('EQSANS_86217', 508339, 1300, 1300, 14122, 9717, 59719, True),
          ('EQSANS_92353', 262291, 4000, 4000, 14122, 11410, 61412, True),
          ('EQSANS_85550', 270022, 5000, 4998, 14122, 12036, 62040, True),
          ('EQSANS_101595', 289989, 1300, 1300, 14122, 7773, 24441, False),
          ('EQSANS_88565', 19362, 4000, 4000, 14122, 45615, 62281, False),
          ('EQSANS_88901', 340431, 8000, 7989, 14122, 67332, 83999, False))
run_sets = [{k: v for k, v in zip(keys, value)} for value in values]


@pytest.fixture(scope='module', params=run_sets)
@namedtuplefy
def rs(request):
    run_set = request.param
    ws = eqsans.load_events(run_set['run'], output_workspace=uwn())
    return {**run_set, **dict(ws=ws)}


class TestLoadEvents(object):

    def test_geometry(self, rs):
        ws = rs.ws
        assert isinstance(ws, EventWorkspace)
        assert ws.getNumberEvents() == rs.num_events
        # assert distance of detector1 same as that in detectorZ of the logs
        instrument = ws.getInstrument()
        det = instrument.getComponentByName(eqsans.detector_name)
        d1 = det.getDistance(instrument.getSample())
        assert rs.nominal_sdd == pytest.approx(d1 * 1000, abs=1)
        # Check logs
        sl = SampleLogs(ws)
        assert rs.ssd == approx(sl.single_value('source-sample-distance'),
                                abs=1)
        assert rs.sdd == approx(sl.single_value('sample-detector-distance'),
                                abs=1)

    def test_tofs(self, rs):
        ws = rs.ws
        assert ws.getTofMin() == pytest.approx(rs.min_tof, abs=1)
        assert ws.getTofMax() == pytest.approx(rs.max_tof, abs=1)
        assert bool(SampleLogs(ws).is_frame_skipping.value) == rs.skip_frame


def test_prepared_data(eqsans_f):
    """
        This is Section 3 of the requirements document.
        It should just be a convenience function that brings together
        loading, moving the detector, normalizing, binning in wavelength,
        and subtracting dark current.
    """
    ws = eqsans.prepare_data(eqsans_f['data'])
    assert isinstance(ws, EventWorkspace)


if __name__ == '__main__':
    pytest.main()
