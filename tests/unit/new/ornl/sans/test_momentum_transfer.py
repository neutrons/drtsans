import numpy as np
from numpy import linalg
import pytest
from mantid.simpleapi import LoadEmptyInstrument, Rebin
from ornl.sans.momentum_transfer import calculate_momentum_transfer


# This implements Issue #210
# dev - Wenduo Zhou <wzz@ornl.gov>
@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 4, 'Ny': 4,
                           'dx': 0.006, 'dy': 0.004, 'zc': 1.25,  # TODO - it is best to use close to real L1 and L2
                           'l1': -5.}],
                         indirect=True)
def test_q_resolution_per_pixel(generic_IDF):
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

    q, qx, qy, qz = calculate_momentum_transfer(ws)
    verify_q(ws, q, qx, qy)

    # Test for multiple bins (wave lengths)
    ws = Rebin(InputWorkspace=ws, Params='2.5, 0.5, 4.5')

    q, qx, qy, qz = calculate_momentum_transfer(ws)
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
    # Obtain wave length
    wavelength_bin_boundary_matrix = ws.extractX()
    wavelength_bin_center_matrix = (wavelength_bin_boundary_matrix[:, 1:] +
                                    wavelength_bin_boundary_matrix[:, :-1]) / 2.0

    # Get instrument information
    spec_info = ws.spectrumInfo()
    num_spec = spec_info.size()

    # Arrays
    ratio_x_vec = np.ndarray(shape=(num_spec, 1), dtype='float')
    ratio_y_vec = np.ndarray(shape=(num_spec, 1), dtype='float')
    pixel_2theta_vec = np.ndarray(shape=(num_spec, 1), dtype='float')

    # sample and moderator information: get K_i
    sample_pos = ws.getInstrument().getSample().getPos()
    source_pos = ws.getInstrument().getSource().getPos()
    k_in = sample_pos - source_pos
    k_in /= linalg.norm(k_in)

    for iws in range(num_spec):
        # get 2theta
        pixel_2theta_vec[iws, 0] = spec_info.twoTheta(iws)  # unit = radians

        # calculate Qx, Qy, Qz
        det_i_pos = ws.getDetector(iws).getPos()
        k_out = det_i_pos - sample_pos
        k_out /= linalg.norm(k_out)

        unit_q = k_out - k_in
        ratio_x_vec[iws] = unit_q[0] / linalg.norm(unit_q)
        ratio_y_vec[iws] = unit_q[1] / linalg.norm(unit_q)
    # END

    # Calculate Q = 4 pi sin(theta)/lambda
    prove_q_matrix = 4.0 * np.pi * np.sin(0.5 * pixel_2theta_vec) / wavelength_bin_center_matrix

    assert np.allclose(q_matrix, prove_q_matrix, rtol=1.E-10)
    assert np.allclose(q_matrix * ratio_x_vec, qx_matrix, rtol=1.E-10)
    assert np.allclose(q_matrix * ratio_y_vec, qy_matrix, rtol=1.E-10)

    return
