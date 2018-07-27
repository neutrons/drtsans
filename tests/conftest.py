from __future__ import (absolute_import, division, print_function)

import sys
import os
import pytest
from os.path import join as pjoin

# Resolve the path to the "external data"
this_module_path = sys.modules[__name__].__file__
parent_dir = pjoin(os.path.dirname(this_module_path), os.pardir)
data_dir = pjoin(parent_dir, 'data')

import mantid.simpleapi as mtds

@pytest.fixture(scope='session')
def eqsans_f():
    return dict(data=pjoin(data_dir, 'eqsans', 'EQSANS_68183_event.nxs'),
                beamcenter=pjoin(data_dir, 'eqsans', 'EQSANS_68200_event.nxs'),
                darkcurrent=pjoin(data_dir, 'eqsans', 'EQSANS_68200_event.nxs')
                )

@pytest.fixture(scope='session')
def eqsans_w(eqsans_f):
    """Load EQSANS files into workspaces"""
    return {k: mtds.LoadEventNexus(v, OutputWorkspace=k)
            for (k, v) in eqsans_f.items()}
