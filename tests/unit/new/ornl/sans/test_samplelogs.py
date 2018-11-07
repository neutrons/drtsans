from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn
from mantid.simpleapi import Load
from mantid.api import Run

from ornl.sans.samplelogs import SampleLogs


class TestSampleLogs(object):

    def test_init(self, refd):
        test_file = pjn(refd.new.sans,
                        'test_samplelogs',
                        'EQSANS_92353_no_events.nxs')
        w = Load(test_file, OutputWorkspace='test_init_w')
        r = w.getRun()
        for other in [test_file, w, r]:
            sl = SampleLogs(other)
            assert isinstance(sl._run, Run)

    def test_getattr(self, refd):
        test_file = pjn(refd.new.sans,
                        'test_samplelogs',
                        'EQSANS_92353_no_events.nxs')
        sl = SampleLogs(test_file)
        assert_almost_equal(sl.Phase1.value.mean(), 22444, decimal=0)


if __name__ == '__main__':
    pytest.main()
