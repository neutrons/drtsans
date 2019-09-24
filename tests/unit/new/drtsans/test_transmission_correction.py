from drtsans.transmission import apply_transmission_correction
import numpy as np
import pytest

# Make an example array of intensity counts
I_sam = np.array([[60,  65,  70,  75,  80],
                  [75,  80,  85,  90,  95],
                  [90,  95, 100, 105, 110],
                  [105, 110, 115, 120, 125],
                  [120, 125, 130, 135, 140]], dtype=np.float)

# Poisson uncertainties on counts
I_sam_err = np.sqrt(I_sam)

# Make an example array of transmission coefficients
transmission = np.array([[0.6, 0.7, 0.7, 0.7, 0.6],
                         [0.7, 0.7, 0.8, 0.7, 0.7],
                         [0.7, 0.8, 0.9, 0.8, 0.7],
                         [0.7, 0.7, 0.8, 0.7, 0.7],
                         [0.6, 0.7, 0.7, 0.7, 0.6]])  # the first entry in this list for some reason is taken as the constant transmission value

# Give the transmission coefficients 5% uncertainties
transmission_err = transmission * 0.05


# Make a mantid workspace for the transmission coefficients
@pytest.mark.parametrize('workspace_with_instrument',
                         [{'Nx': I_sam.shape[0], 'Ny': I_sam.shape[1]}],
                         indirect=True)
def test_transmission_correction(workspace_with_instrument):
    '''
    Test the calculation of the detector transmission correction and error propagation
    given in Eq. 7.7 and 7.8 in the master document
    dev - Joe Osborn <osbornjd@ornl.gov>
    '''

    # Create the workspaces
    I_sam_wksp = workspace_with_instrument(axis_values=[5.95, 6.075],
                                           intensities=I_sam,
                                           uncertainties=I_sam_err)

    transmission_wksp = workspace_with_instrument(axis_values=[5.95, 6.075],
                                                  intensities=transmission,
                                                  uncertainties=transmission_err)

    # Check that they were created successfully
    assert I_sam_wksp.extractY().sum() == pytest.approx(np.sum(I_sam))
    assert transmission_wksp.extractY().sum() == pytest.approx(np.sum(transmission))
    assert I_sam_wksp.extractE().sum() == pytest.approx(np.sum(I_sam_err))
    assert transmission_wksp.extractE().sum() == pytest.approx(np.sum(transmission_err))

    # Calculate the result by hand to compare to drtsans
    # Calculate transmission corrected intensity
    I_tsam = I_sam / transmission

    # Propagate uncertainties to transmission corrected intensity
    I_tsam_err = I_tsam * np.sqrt((I_sam_err/I_sam)**2 + (transmission_err/transmission)**2)

    # Calculate the result with drtsans framework
    result = apply_transmission_correction(I_sam_wksp, transmission_wksp)

    # Compare the result from drtsans and the result calcualted by hand
    np.testing.assert_equal(result.extractY().reshape(I_tsam.shape), I_tsam)
    np.testing.assert_equal(result.extractE().reshape(I_tsam.shape), I_tsam_err)


if __name__ == '__main__':
    pytest.main()
