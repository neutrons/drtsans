import pytest
import numpy as np
from drtsans.transmission import apply_transmission_correction


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='GPSANS', Nx=9, Ny=1, dx=0.1, dy=0.01, zc=0.3)], indirect=True)
def test_calculate_theta_dependent_transmission_single_value(workspace_with_instrument):
    r"""
    This implements Issue #176, addressing master document section 7
    dev - Jose Borreguero <borreguerojm@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>

    The instrument is made up of 9 tubes, each tube containing only one square pixel of side 10cm.
    The instrument is positioned 30cm away from the sample. These specifications result in two_theta
    angles matching those of the test assembled by the instrument team.
    """
    transmission_zero_angle, transmission_zero_angle_error = 0.9, 0.1

    # Input sample workspace
    wavelength_bin = [3.0, 3.2]  # some arbitrary wavelength bin
    sample_counts = np.ones(9).reshape((9, 1))  # each pixel has an intensity of one
    sample_errors = np.zeros(9).reshape((9, 1))  # only interested in errors arising in the transmission calculation
    sample_workspace = workspace_with_instrument(axis_values=wavelength_bin, intensities=sample_counts,
                                                 uncertainties=sample_errors, view='array')

    # Apply the transmission
    sample_workspace = apply_transmission_correction(sample_workspace, trans_value=transmission_zero_angle,
                                                     trans_error=transmission_zero_angle_error, theta_dependent=True)

    # Manually obtain the transmission for each detector
    pixel_ids = sample_workspace.getNumberHistograms()
    spectrum_info = sample_workspace.spectrumInfo()
    two_thetas = np.array([spectrum_info.twoTheta(pixel_id) for pixel_id in range(pixel_ids)])
    # pixel
    exponents = 0.5 / np.cos(two_thetas) + 0.5
    manual_transmissions = np.power(transmission_zero_angle, exponents)
    manual_transmissions_errors = exponents * manual_transmissions * transmission_zero_angle_error / \
        transmission_zero_angle  # errors in two_theta are not included

    # Compare transmissions from the sample workspace with those calculated manually
    transmissions = 1.0 / sample_workspace  # recall sample_workspace contains corrected intensities
    assert transmissions.extractY().flatten() == pytest.approx(manual_transmissions, abs=0.001)
    assert transmissions.extractE().flatten() == pytest.approx(manual_transmissions_errors, abs=0.001)
