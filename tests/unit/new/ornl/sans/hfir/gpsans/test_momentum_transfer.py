#!/usr/bin/env python
from __future__ import print_function

import tempfile

import numpy as np

from mantid import mtd
from mantid.simpleapi import Load, LoadHFIRSANS, Rebin2D
from ornl.sans.momentum_transfer import (bin_into_q1d, bin_into_q2d,
                                         bin_wedge_into_q1d)
from reduction_workflow.command_interface import AppendDataFile, Reduce
from reduction_workflow.instruments.sans.hfir_command_interface import (
    GPSANS, AzimuthalAverage, IQxQy, OutputPath, SetBeamCenter)


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
    '''

    '''

    GPSANS()
    SetBeamCenter(92, 123)
    AppendDataFile(gpsans_f['sample_scattering_2'], workspace='ws_legacy')
    # AppendDataFile(
    #     '/home/rhf/git/sans-rewrite/data/new/ornl/sans/hfir/gpsans/\
    # CG2_exp325_scan0007_0001.xml', workspace='ws_legacy')
    AzimuthalAverage(binning="0.01,0.001,0.11")
    OutputPath(tempfile.gettempdir())
    IQxQy(nbins=110, log_binning=False)
    Reduce()

    ws_legacy = mtd['ws_legacy']
    ws_legacy_iq = mtd['ws_legacy_iq']
    ws_legacy_iqxqy = mtd['ws_legacy_iqxy']

    wss_name_ws = bin_into_q2d(ws_legacy, out_ws_prefix='ws_new')
    ws_new_iqxqy, ws_new_dqx, ws_new_dqy = [ws[1] for ws in wss_name_ws]

    ws_new_iqxqy = Rebin2D(ws_new_iqxqy, Axis1Binning='-0.11,0.001,0.11',
                           Axis2Binning='-0.11,0.001,0.11')
    ws_new_dqx = Rebin2D(ws_new_dqx, Axis1Binning='-0.11,0.001,0.11',
                         Axis2Binning='-0.11,0.001,0.11')
    ws_new_dqy = Rebin2D(ws_new_dqy, Axis1Binning='-0.11,0.001,0.11',
                         Axis2Binning='-0.11,0.001,0.11')

    legacy_iqxqy = ws_legacy_iqxqy.extractY()
    new_iqxqy = ws_new_iqxqy.extractY()
    # this is failling
    assert np.allclose(legacy_iqxqy, new_iqxqy, atol=1e-01)

    _, ws_new_iq = bin_into_q1d(ws_new_iqxqy, ws_new_dqx, ws_new_dqy,
                                bins=np.linspace(0.01, 0.11, 101),
                                out_ws_prefix='ws_new')

    legacy_iq = ws_legacy_iq.extractY()
    new_iq = ws_new_iq.extractY()

    legacy_iq = np.nan_to_num(legacy_iq[0])
    new_iq = np.nan_to_num(new_iq[0])

    assert np.allclose(legacy_iq, new_iq, atol=1e-01)
