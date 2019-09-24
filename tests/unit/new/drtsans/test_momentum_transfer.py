import numpy as np
import pytest
from mantid.simpleapi import LoadEmptyInstrument, Rebin
from drtsans.momentum_transfer import calculate_momentum_transfer


# This implements Issue #210: calculate general algorithm to calculate Q, Qx and Qy
# dev - Wenduo Zhou <wzz@ornl.gov>
@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4,
                           'dx': 0.006, 'dy': 0.004, 'zc': 1.25,
                           'l1': -5.}],
                         indirect=True)
def test_calculate_momentum_transfer(generic_IDF):
    """
    Create a generic (SANS) instrument and test for main()
    :param generic_IDF:
    :return:
    """
    # Generate a generic SANS instrument with a pixel of
    # the size and position specified in
    # sans-backend/documents/Master_document_022219.pdf
    with open(r'/tmp/GenericSANS_Definition.xml', 'w') as tmp:
        tmp.write(generic_IDF)
        tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS',
                             OutputWorkspace='test_uncertainty')
    ws.getAxis(0).setUnit('Wavelength')

    # Test for single bin
    # Range of TOF
    wave_length_range = np.array([2.5, 4.5])  # A
    intensity = np.array([20.]*16) + np.random.randn(16,)
    init_delta_intensity = np.random.randn(intensity.shape[0], )

    # assume that the TOF is already frame corrected
    for i in range(16):
        ws.dataX(i)[:] = wave_length_range  # microseconds
        ws.dataY(i)[0] = intensity[i]
        ws.dataE(i)[0] = init_delta_intensity[i]

    # TODO FIXME - Need to review because Qx and Qy are redefined
    q, qx, qy, two_theta, s2p = calculate_momentum_transfer(ws)
    verify_q(ws, q, qx, qy)

    # Test for multiple bins (wave lengths)
    ws = Rebin(InputWorkspace=ws, Params='2.5, 0.5, 4.5')

    q, qx, qy, qz, s2p = calculate_momentum_transfer(ws)
    verify_q(ws, q, qx, qy)

    return


def verify_q(ws, q_matrix, qx_matrix, qy_matrix):
    """
    Verify Q, Qx and Qy values calculated against Q values calculated from other algorithms that
    calculates Q with 4 pi sin(theta) / lambda, Qx and Qy from Q
    :param ws: workspace instance
    :param q_matrix: N x M array (N is number of histogram, M is number of Q values)
    :param qx_matrix: N x M array (N is number of histogram, M is number of Q values)
    :param qy_matrix: N x M array (N is number of histogram, M is number of Q values)
    :return:
    """
    # Check dimension
    num_spec = ws.getNumberHistograms()
    num_bins = ws.readY(0)[0]

    assert q_matrix.shape == (num_spec, num_bins)
    assert qx_matrix.shape == (num_spec, num_bins)
    assert qy_matrix.shape == (num_spec, num_bins)

    # Test value
    assert np.allclose(q_matrix**2, qx_matrix**2 + qy_matrix**2, rtol=1.E-10)

    return
