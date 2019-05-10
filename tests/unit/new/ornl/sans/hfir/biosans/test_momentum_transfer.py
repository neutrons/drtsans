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

    wss_name_ws = bin_into_q2d(ws, "detector1", out_ws_prefix="main")
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    wss_name_ws = bin_into_q2d(ws, "wing_detector", out_ws_prefix="wing")
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 160)
    assert ws_iqxqy.extractX().shape == (256, 161)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)


def bin_into_q2d_parallel(parameters):
    ws_name, component_name, out_ws_prefix = parameters
    ws = mtd[ws_name] # need to pass the name. ws is shared between 2 processes?
    workspaces = bin_into_q2d(ws, component_name, out_ws_prefix)
    return workspaces


def bin_into_q1d_parallel(parameters):
    ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix = parameters
    iq_ws = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix=out_ws_prefix)
    return iq_ws


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
        results = pool.map(bin_into_q2d_parallel, parameters)

    list_of_names_and_wss = [ws for ws_components in results
                             for ws in ws_components]
    # Put the wss results in the AnalysisDataService
    [AnalysisDataService.add(ws[0], ws[1]) for ws in list_of_names_and_wss]
    # see if those wss are a sublist in AnalysisDataService
    assert all(elem in AnalysisDataService.getObjectNames() for elem in [
        ws[0] for ws in list_of_names_and_wss])

    # list with (ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix)
    parameters = [[ws[1] for ws in ws_components] + [out_ws_prefix]
                  for ws_components, out_ws_prefix
                  in zip(results, ['main', 'wing'])]

    with multiprocessing.Pool(processes=2) as pool:
        results = pool.map(bin_into_q1d_parallel, parameters)

    [AnalysisDataService.add(ws[0], ws[1]) for ws in results]
    assert all(elem in AnalysisDataService.getObjectNames() for elem in [
        ws[0] for ws in results])
