from drtsans.transmission import apply_transmission_correction
import numpy as np
import matplotlib.pyplot as plt
import pytest

I_sam = np.array([[ 60,  65,  70,  75,  80],
                  [ 75,  80,  85,  90,  95],
                  [ 90,  95, 100, 105, 110],
                  [105, 110, 115, 120, 125],
                  [120, 125, 130, 135, 140]])

transmission = np.array( [[0.6, 0.7, 0.7, 0.7, 0.6],
                          [0.7, 0.7, 0.8, 0.7, 0.7],
                          [0.7, 0.8, 0.9, 0.8, 0.7],
                          [0.7, 0.7, 0.8, 0.7, 0.7],
                          [0.6, 0.7, 0.7, 0.7, 0.6]] )

# Give the transmission coefficients 5% uncertainties
transmission_err = np.array(transmission * 0.05)
    
pixel_size = 0.005 # meters

# Make a mantid workspace for the intensities
@pytest.mark.parametrize('generic_workspace',
                         [{'name': 'I_sam',
                           'dx': pixel_size, 'dy': pixel_size,
                           'yc': pixel_size / 2.,
                           'axis_values': [5.95,6.075],
                            'intensities': I_sam}],
                         indirect = True)

# Make a mantid workspace for the transmission coefficients
@pytest.mark.parametrize('generic_workspace',
                          [{'name': 'transmission',
                            'intensities' : transmission,
                            'uncertainties' : transmission_err
                            }],
                          indirect=True)

# I can get a single workspace through to the function, but when I try to use
# two workspaces I run into issues
def test_transmission_correction(generic_workspace, transmission_workspace):
    '''
    Test the calculation of the detector transmission correction and error propagation
    given in Eq. 7.7 and 7.8 in the master document
    dev - Joe Osborn <osbornjd@ornl.gov>
    '''
    
    I_sam_wksp = generic_workspace
    transmission_wksp = transmission_workspace;
    
    # Calculate the result with drtsans framework
    result = apply_transmission_correction(I_sam_wksp, transmission_wksp, 1, 0.05)

    # Calculate the result by hand to compare to drtsans
    # Poisson uncertainties on counts
    I_sam_err = np.sqrt(I_sam)

    # Calculate transmission corrected intensity
    I_tsam = I_sam / transmission

    # Propagate uncertainties to transmission corrected intensity
    I_tsam_err = I_sam / transmission * np.sqrt( (I_sam_err/I_sam)**2 + (transmission_err/transmission)**2 )

    # Compare the result from drtsans and the result calcualted by hand
    assert result.extractY()[0][0] == I_tsam[0][0]
    assert result.extractE()[0][0] == I_tsam_err[0][0]
 
if __name__ == '__main__':
    pytest.main()
  
