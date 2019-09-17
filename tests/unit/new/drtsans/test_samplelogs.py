from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn
from mantid.simpleapi import Load
from mantid.api import Run

from drtsans.samplelogs import SampleLogs


@pytest.mark.offline
class TestSampleLogs(object):
    def test_init(self, reference_dir):
        test_file = pjn(reference_dir.new.sans,
                        'test_samplelogs',
                        'EQSANS_92353_no_events.nxs')
        w = Load(test_file, OutputWorkspace='test_init_w')
        r = w.getRun()
        for other in [test_file, w, r]:
            sl = SampleLogs(other)
            assert isinstance(sl._run, Run)

    def test_getattr(self, reference_dir):
        test_file = pjn(reference_dir.new.sans,
                        'test_samplelogs',
                        'EQSANS_92353_no_events.nxs')
        sl = SampleLogs(test_file)
        assert_almost_equal(sl.Phase1.value.mean(), 22444, decimal=0)

    def test_insert(self, reference_dir):
        test_file = pjn(reference_dir.new.sans,
                        'test_samplelogs',
                        'EQSANS_92353_no_events.nxs')

        sl = SampleLogs(test_file)
        sl.insert('string_log', 'log value')
        assert sl.string_log.value, 'log value'
        assert not sl.string_log.units

        units = 'super awesome units'
        sl.insert('int_log', 42, units)
        assert sl.int_log.value == 42
        assert sl.int_log.units == units

        euler = 2.7182818284590452353602874713527
        units = 'even more awesome units'
        sl.insert('float_log', euler, units)
        assert sl.float_log.units == units
        assert sl.float_log.value == euler

        values = list(range(1, 9))
        units = 'most awesomest units ever'
        sl.insert('array_log', values, units)
        assert sl.array_log.units == units
        # this seems like the wrong value from mantid
        assert sl.array_log.value == values[0]


if __name__ == '__main__':
    pytest.main()
