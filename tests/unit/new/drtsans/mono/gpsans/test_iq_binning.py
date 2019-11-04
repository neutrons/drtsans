#!/usr/bin/env python
import numpy as np
import scipy

# from drtsans.iq import bin_iq_into_linear_q2d, bin_iq_into_linear_q1d, BinningMethod
from mantid.simpleapi import LoadHFIRSANS

# This test implements issue #??? to verify methods to bin I(Q) for GP-SANS
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Lisa???


def skip_test_momentum_tranfer_with_annular_1d_binning(gpsans_f):
    '''
    Test Iq, Iqxqy with the ends of the detector tubes masked
    '''
    filename = gpsans_f['sample_scattering_2']
    ws = LoadHFIRSANS(Filename=filename)

    mt = calculate_momentum_transfer(ws)  # noqa: F821

    _, ws_iq = mt.bin_into_q1d()
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    q_min = 0.05
    q_max = 0.2
    _, ws_annular_iq = mt.bin_annular_into_q1d(q_min=q_min,
                                               q_max=q_max,
                                               bins=20)
    assert ws_annular_iq.extractY().shape == (1, 20)
    assert ws_annular_iq.extractX().shape == (1, 21)

    # we are comparing the I of both curves. Since the binning and number of
    # bins is different, we need to interpolate the values
    iq_q = ws_iq.extractX().ravel()
    iq_q = (iq_q[1:] + iq_q[:-1]) / 2.0  # Bin centres
    iq_q_subset = iq_q[(iq_q >= q_min) & (iq_q <= q_max)]
    iq_i = ws_iq.extractY().ravel()
    iq_i_subset = iq_i[(iq_q >= q_min) & (iq_q <= q_max)]

    iq_annular_q = ws_annular_iq.extractX().ravel()
    iq_annular_q = (iq_annular_q[1:] + iq_annular_q[:-1]) / 2.0  # Bin centres
    iq_annular_i = ws_annular_iq.extractY().ravel()

    f = scipy.interpolate.InterpolatedUnivariateSpline(iq_annular_q,
                                                       iq_annular_i)

    assert np.allclose(iq_i_subset, f(iq_q_subset), rtol=1)
