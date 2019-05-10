#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadHFIRSANS
from ornl.settings import unique_workspace_name
from ornl.sans.momentum_transfer import bin_into_q2d, bin_into_q1d
from mantid import mtd

def test_momentum_tranfer(biosans_sensitivity_dataset):

    ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'],
        OutputWorkspace=unique_workspace_name())

    qxqy_wss_grouped = bin_into_q2d(ws, "detector1", out_ws_prefix="main")
    assert len(qxqy_wss_grouped) == 3

    ws_iqxqy, ws_dqx, ws_dqy = qxqy_wss_grouped
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    qxqy_wss_grouped = bin_into_q2d(ws, "wing_detector", out_ws_prefix="wing")
    assert len(qxqy_wss_grouped) == 3

    ws_iqxqy, ws_dqx, ws_dqy = qxqy_wss_grouped
    assert ws_iqxqy.extractY().shape == (256, 160)
    assert ws_iqxqy.extractX().shape == (256, 161)

    ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)



def parallel_f(parameters):
    ws_name, component_name, out_ws_prefix = parameters
    ws = mtd[ws_name]
    qxqy_wss_grouped = bin_into_q2d(ws, component_name, out_ws_prefix)
    return [(ws.name, ws) for ws in qxqy_wss_grouped]

def test_momentum_tranfer_parallel(biosans_sensitivity_dataset):

    import multiprocessing
    from mantid.simpleapi import AnalysisDataService

    input_ws_name = unique_workspace_name()
    LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'],
        OutputWorkspace=input_ws_name)
    
    ws_names = [input_ws_name, input_ws_name]
    components = ['detector1', 'wing_detector']
    prefixes = ['main', 'wing']
    parameters = [
        (ws_name, component, prefix) for ws_name, component, prefix in
        zip(ws_names, components, prefixes)
    ]
    with multiprocessing.Pool(processes=2) as pool:
        results = pool.map(parallel_f, parameters)

    [AnalysisDataService.add(name, ws) for name, ws in zip([
        ('main_iqxqy', 'main_dqx', 'main_dqy'),
        ('wing_iqxqy', 'wing_dqx', 'wing_dqy')], results)]

