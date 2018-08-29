from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal

from ornl.sans.sns.eqsans import tofcorrect


def test_frame_width(porasil_slice1m):
    s = porasil_slice1m.w['s']  # sample workspace
    fw = tofcorrect.frame_width(s)
    assert_almost_equal(fw, 16666, 0)


def test_frame_skipping(porasil_slice1m, frame_skipper):
    s = porasil_slice1m.w['s']  # sample workspace
    assert tofcorrect.frame_skipping(s) is False
    s = frame_skipper.w['s']  # sample workspace
    assert tofcorrect.frame_skipping(s) is True


if __name__ == '__main__':
    pytest.main()
