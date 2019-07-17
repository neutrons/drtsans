#!/usr/bin/env python
from __future__ import print_function

import tempfile

import numpy as np

import mantid
from mantid import mtd
from mantid.simpleapi import LoadHFIRSANS, MaskBTP, CloneWorkspace
from ornl.sans.momentum_transfer import MomentumTransfer
from reduction_workflow.command_interface import AppendDataFile, Reduce
from reduction_workflow.instruments.sans.hfir_command_interface import (
    GPSANS, AzimuthalAverage, IQxQy, OutputPath, SetBeamCenter)


def test_momentum_tranfer_wedge_anisotropic(gpsans_f):
    '''
    Tests the basic for momentum transfer:
    - Iq, Iqxqy: Shape of the output workspaces
    - Wedge location: makes sure the feature has a higher count
    '''

    ws = LoadHFIRSANS(
        Filename=gpsans_f['anisotropic'],
        OutputWorkspace='aniso_raw')

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

    _, ws_iq_feature = mt.bin_wedge_into_q1d(phi_0=55, phi_aperture=30)
    ws_iq_feature_i = ws_iq_feature.extractY()
    assert ws_iq_feature_i.shape == (1, 100)
    assert ws_iq_feature.extractX().shape == (1, 101)
    ws_iq_feature = CloneWorkspace(ws_iq_feature)

    _, ws_iq_non_feature = mt.bin_wedge_into_q1d(phi_0=55+90, phi_aperture=30)
    ws_iq_non_feature_i = ws_iq_non_feature.extractY()
    assert ws_iq_non_feature_i.shape == (1, 100)
    assert ws_iq_non_feature.extractX().shape == (1, 101)
    ws_iq_non_feature = CloneWorkspace(ws_iq_non_feature)

    # Wedge with feature has more counts than that without it
    assert (ws_iq_feature_i.sum() - ws_iq_non_feature_i.sum()) > 1000


def test_momentum_tranfer_cross_check(gpsans_f):
    '''
    Tests Iq on the legacy vs new reduction
    '''

    GPSANS()
    SetBeamCenter(92, 123)
    AppendDataFile(
        gpsans_f['sample_scattering_2'],
        workspace='ws_legacy')
    AzimuthalAverage(binning="0.01,0.001,0.11")
    OutputPath(tempfile.gettempdir())
    IQxQy(nbins=110, log_binning=False)
    Reduce()

    ws_legacy = mtd['ws_legacy']
    ws_legacy_iq = mtd['ws_legacy_iq']

    mt = MomentumTransfer(ws_legacy, out_ws_prefix='ws_new')

    _, ws_new_iq = mt.bin_into_q1d(bins=np.linspace(0.01, 0.11, 101))

    legacy_iq = ws_legacy_iq.extractY()
    new_iq = ws_new_iq.extractY()

    legacy_iq = np.nan_to_num(legacy_iq[0])
    new_iq = np.nan_to_num(new_iq[0])

    assert np.allclose(legacy_iq, new_iq)


def test_momentum_tranfer_with_and_without_mask(gpsans_f):
    '''
    Test Iq, Iqxqy with the ends of the detector tubes masked
    '''
    filename = gpsans_f['sample_scattering_2']
    ws_mask = LoadHFIRSANS(Filename=filename, OutputWorkspace='mask_raw')
    ws_no_mask = LoadHFIRSANS(Filename=filename, OutputWorkspace='no_mask_raw')

    # Let's mask the detector ends: 20 + 20
    MaskBTP(ws_mask, Components='detector1', Pixel='0-19,236-255')

    mt_mask = MomentumTransfer(ws_mask, out_ws_prefix='mask')
    mt_no_mask = MomentumTransfer(ws_no_mask, out_ws_prefix='no_mask')

    _, ws_mask = mt_mask.bin_into_q2d()
    _, ws_no_mask = mt_no_mask.bin_into_q2d()

    # ws_no_mask has more Q range because it was not masked
    for i in range(10, 19):
        assert ws_mask.readY(i).any() != ws_no_mask.readY(i).any()

    # Let's make sure the I(qx,qy) has masked the values as the detector
    i = mt_mask.i.filled()
    nans = np.isnan(i)
    assert np.any(nans)
    assert np.count_nonzero(nans) == 192 * (20 + 20)  # 192 tubes

    _, ws_mask_iq = mt_mask.bin_into_q1d()
    assert ws_mask_iq.extractY().shape == (1, 100)
    assert ws_mask_iq.extractX().shape == (1, 101)

    _, ws_non_mask_iq = mt_no_mask.bin_into_q1d()
    assert ws_non_mask_iq.extractY().shape == (1, 100)
    assert ws_non_mask_iq.extractX().shape == (1, 101)
    # Q range in non mask is higher
    assert ws_non_mask_iq.extractX().max() > ws_mask_iq.extractX().max()
