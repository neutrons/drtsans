#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

import pytest
from mantid.simpleapi import Load

from ornl.settings import amend_config, unique_workspace_name
from ornl.sans.sns.eqsans.beam_finder import direct_beam_center


def test_direct_beam_center():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        w = Load(Filename='EQSANS_92160',
                 OutputWorkspace=unique_workspace_name())
    assert direct_beam_center(w) == pytest.approx((0.025, 0.013), abs=1e-3)


if __name__ == '__main__':
    pytest.main()
