#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadHFIRSANS, Load
from ornl.sans.momentum_transfer import (
    bin_into_q1d, bin_into_q2d, bin_wedge_into_q1d)


def test_momentum_tranfer_anisotropic(gpsans_f):

    ws = LoadHFIRSANS(
        Filename=gpsans_f['anisotropic'],
        OutputWorkspace='aniso_raw')

    wss_name_ws = bin_into_q2d(ws, out_ws_prefix='aniso')
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy, out_ws_prefix='aniso')
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    _, ws_iq_feature = bin_wedge_into_q1d(ws_iqxqy, ws_dqx, ws_dqy,
                                          phi_0=55, phi_aperture=30,
                                          out_ws_prefix='ws_feature')
    ws_iq_feature_i = ws_iq_feature.extractY()
    assert ws_iq_feature_i.shape == (1, 100)
    assert ws_iq_feature.extractX().shape == (1, 101)

    _, ws_iq_non_feature = bin_wedge_into_q1d(ws_iqxqy, ws_dqx, ws_dqy,
                                              phi_0=55+90, phi_aperture=30,
                                              out_ws_prefix='ws_non_feature')
    ws_iq_non_feature_i = ws_iq_non_feature.extractY()
    assert ws_iq_non_feature_i.shape == (1, 100)
    assert ws_iq_non_feature.extractX().shape == (1, 101)

    # Wedge with feature has more counts than that without it
    assert (ws_iq_feature_i.sum() - ws_iq_non_feature_i.sum()) > 1000


def test_momentum_tranfer_cross_check(gpsans_f):

    ws = Load(Filename=gpsans_f['ws'], OutputWorkspace='ws')
    ws_iq = Load(Filename=gpsans_f['ws_iq'], OutputWorkspace='ws_iq')
    ws_iqxqy = Load(Filename=gpsans_f['ws_iqxqy'], OutputWorkspace='ws_iqxqy')

    # This fails !!
    # wss_name_ws = bin_into_q2d(ws, out_ws_prefix='ws')
    # ws_iqxqy_reduced = wss_name_ws[0][1]
