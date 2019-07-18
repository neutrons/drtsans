#!/usr/bin/env python
from __future__ import print_function

from multiprocessing import Pool

import pytest

import mantid
from mantid import mtd
from mantid.kernel import logger
from mantid.simpleapi import CloneWorkspace, LoadHFIRSANS
from ornl.sans.momentum_transfer import MomentumTransfer
from ornl.settings import unique_workspace_name


def test_momentum_tranfer_serial(biosans_f):

    ws = LoadHFIRSANS(
        # Filename='/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/'\
        #      'biosans/BioSANS_exp440_scan0022_0006.xml',
        Filename=biosans_f['anisotropic'],
        OutputWorkspace=unique_workspace_name())

    mt = MomentumTransfer(ws)
    assert mt.qx.shape == mt.qy.shape == mt.dqx.shape == mt.dqy.shape == \
        mt.i.shape == mt.i_sigma.shape == (256*192, )

    table_iq = mt.q2d()
    assert type(table_iq) == mantid.dataobjects.TableWorkspace

    _, ws = mt.bin_into_q2d()
    assert ws.extractY().shape == (256, 192)
    assert ws.extractX().shape == (256, 193)

    _, ws = mt.bin_into_q1d()
    assert ws.extractY().shape == (1, 100)
    assert ws.extractX().shape == (1, 101)

    _, ws_iq_feature = mt.bin_wedge_into_q1d(phi_0=0, phi_aperture=30)
    ws_iq_feature_i = ws_iq_feature.extractY()
    assert ws_iq_feature_i.shape == (1, 100)
    assert ws_iq_feature.extractX().shape == (1, 101)
    ws_iq_feature = CloneWorkspace(ws_iq_feature)

    _, ws_iq_non_feature = mt.bin_wedge_into_q1d(phi_0=0+90, phi_aperture=30)
    ws_iq_non_feature_i = ws_iq_non_feature.extractY()
    assert ws_iq_non_feature_i.shape == (1, 100)
    assert ws_iq_non_feature.extractX().shape == (1, 101)
    ws_iq_non_feature = CloneWorkspace(ws_iq_non_feature)

    # Wedge with feature has more counts than that without it
    assert (ws_iq_feature_i.sum() - ws_iq_non_feature_i.sum()) > 200


def bin_in_parallel(params):

    component_name, out_ws_prefix = params
    ws = mtd["ws_data_raw"]
    mt = MomentumTransfer(ws, component_name=component_name,
                          out_ws_prefix=out_ws_prefix)

    table_iq = mt.q2d()

    ws_q2d_name, ws = mt.bin_into_q2d()

    ws_q1d_name, ws = mt.bin_into_q1d()

    logger.notice("Returning from the parallel calculation: {} {} {}".format(
        table_iq.name(), ws_q2d_name, ws_q1d_name))
    logger.notice("Workspaces names in the parallel function: {}".format(
        mtd.getObjectNames()))

    return table_iq.name(), ws_q2d_name, ws_q1d_name


@pytest.mark.skip(reason="Does not create the WSs in parallel."
                  "Waiting for the mantid team solve the bug")
def test_momentum_tranfer_parallel(biosans_f):
    '''
    For now let's skip this test. it does not run as part of multiple tests
    It runs as a single test though:
    pytest tests/unit/new/ornl/sans/hfir/biosans/test_momentum_transfer.py::\
        test_momentum_tranfer_parallel
    '''

    ws = LoadHFIRSANS(
        # Filename='/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/'
        #     'biosans/BioSANS_exp440_scan0022_0006.xml',
        Filename=biosans_f['anisotropic'],
        OutputWorkspace="ws_data_raw")

    parameters = [('detector1', 'main'), ('wing_detector', 'wing')]
    with Pool(processes=2) as pool:
        return_list = pool.map(bin_in_parallel, parameters)

    logger.notice("After launching return_list: {}".format(return_list))
    logger.notice("After launching we have the WSs: {}".format(
        mtd.getObjectNames()))

    for ret in return_list:
        if ret[0].startswith('wing'):
            y_shape = (256, 160)
            x_shape = (256, 161)
        else:
            y_shape = (256, 192)
            x_shape = (256, 193)

        ws = mtd[ret[1]]
        assert ws.extractY().shape == y_shape
        assert ws.extractX().shape == x_shape

        ws = mtd[ret[2]]
        assert ws.extractY().shape == (1, 100)
        assert ws.extractX().shape == (1, 101)
