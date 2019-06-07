#!/usr/bin/env python
from __future__ import print_function

from mantid import mtd
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.settings import unique_workspace_name

'''
import pytest
import os
os.chdir("/home/rhf/git/sans-rewrite")

pytest.main(["-vs", "/home/rhf/git/sans-rewrite/tests/integration/new/ornl/\
sans/sns/eqsans/test_momentum_transfer.py"])
'''


def bin_into_q2d_parallel(parameters):
    ws_name, component_name, out_ws_prefix = parameters
    ws = mtd[ws_name]  # need to pass the name. ws is shared between processes?
    workspaces = bin_into_q2d(ws, component_name, out_ws_prefix)
    return workspaces


def bin_into_q1d_parallel(parameters):
    ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix = parameters
    iq_ws = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix=out_ws_prefix)
    return iq_ws


# @pytest.mark.skip(reason="It only passes if run as standalone test")1
def test_momentum_tranfer_parallel():
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel

    TODO!
    '''

    import multiprocessing
    import numpy as np
    from mantid.simpleapi import AnalysisDataService, Rebin, ExtractSpectra
    from ornl.sans.sns.eqsans import reduce
    ws = reduce.load_w('EQSANS_68200', unique_workspace_name(suffix="_raw"),
                       low_tof_clip=500, high_tof_clip=2000, dw=0.1)
    ws = Rebin(
        InputWorkspace=ws,
        OutputWorkspace=unique_workspace_name(suffix="_rebin"),
        Params='2.6,0.2,5.6')

    bins = np.arange(2.6, 5.6, 0.2)

    for index, bin_start in enumerate(bins):
        ws_extract = ExtractSpectra(
            InputWorkspace=ws,
            OutputWorkspace=unique_workspace_name(suffix="_{:}".format(index)),
            XMin=bin_start, XMax=bin_start+0.2)
