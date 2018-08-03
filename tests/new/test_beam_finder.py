#!/usr/bin/env python
from __future__ import print_function


import pytest


def test_beam_finder():
    '''
    Test with the new beam finder
    '''

    from ornl.sans.sns.eqsans import beam_finder

    x, y = beam_finder.direct_beam_center()

    assert x == pytest.approx(0.02652545)
    assert y == pytest.approx(0.01804158)
