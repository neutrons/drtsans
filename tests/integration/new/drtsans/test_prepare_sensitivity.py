import numpy as np
import copy
from mantid.simpleapi import MaskDetectors
import pytest
from drtsans.sensitivity import prepare_sensitivity_correction
from numpy.testing import assert_allclose

@pytest.mark.parametrize('workspace_with_instrument',
                          [dict(name='EQSANS', Nx=20, Ny=8)], indirect=True)
def test_prepare_sensitivity(workspace_with_instrument):

    flood_field_measurement = np.array([[65, 68, 66, 75, 71,  68, 66, 70],
                                        [69, 65, 69, 71, 71,  68, 68, 66],
                                        [75, 69, 70, 67, 66,  74, 71, 70],
                                        [66, 71, 71, 70, 71,  68, 66, 65],
                                        [67, 69, 66, 74, 75,  72, 71, 68],
                                        [73, 69, 68, 71, 66,  72, 70, 73],
                                        [67, 6,	 66, 74, 74,  65, 65, 70],
                                        [72, 67, 69, 71, 74,  68, 75, 68],
                                        [71, 75, 72, 73, 73,  69, 69, 67],
                                        [69, 69, 70, 66, 72,  72, 67, 74],
                                        [72, 75, 70, 67, 74,  67, 68, 73],
                                        [69, 66, 67, 68, 68,  75, 71, 73],
                                        [70, 71, 73, 74, 67, 492, 70, 75],
                                        [73, 65, 72, 66, 70,  67, 66, 70],
                                        [69, 68, 71, 68, 70,  72, 67, 70],
                                        [66, 72, 69, 70, 66,  66, 70, 74],
                                        [65, 65, 67, 72, 69,  75, 75, 73],
                                        [65, 72, 72, 75, 67,  73, 75, 72],
                                        [67, 65, 69, 71, 68,  65, 71, 70],
                                        [72, 72, 65, 75, 68,  74, 75, 71]])

    flood_field_measurement_uncertainty = np.sqrt(flood_field_measurement)

    mask = np.ones((20, 8))
    mask[0, :] = np.nan
    mask[8, 3] = np.nan
    mask[9, 3] = np.nan
    mask[19, :] = np.nan

    ffm_with_mask = mask*flood_field_measurement
    ffm_uncertainty_with_mask = mask*flood_field_measurement_uncertainty
    F = np.nanmean(ffm_with_mask)
    n_elements = np.sum(np.logical_not(np.isnan(ffm_with_mask)))
    dF = np.sqrt(np.nansum(np.power(ffm_uncertainty_with_mask, 2)))/n_elements
    II = ffm_with_mask/F
    dI = II * np.sqrt(np.square(ffm_uncertainty_with_mask/ffm_with_mask) + np.square(dF/F))

    interp = np.array([[-5.55e-4,  1.3720e-2, 0.892143],
                       [-6.55e-4,  1.2996e-2, 0.909765],
                       [ -8.9e-5,  4.72e-4, 0.967609],
                       [ 2.96e-4, -5.991e-3, 0.998240],
                       [-6.63e-4,  1.5604e-2, 0.899279],
                       [0.000000, -6.4e-5, 0.969006],
                       [ 4.34e-4, -1.0815e-2, 1.017307],
                       [-5.71e-4,  6.709e-3, 0.980341]])

    interp_uncertainty = np.array([[4.01e-4, 7.882e-3, 0.032903],
                                   [4.99e-4, 9.798e-3, 0.040252],
                                   [3.19e-4, 6.248e-3, 0.025796],
                                   [4.72e-4, 9.217e-3, 0.035956],
                                   [4.10e-4, 8.021e-3, 0.033378],
                                   [5.53e-4, 1.0765e-2, 0.043732],
                                   [4.40e-4, 8.559e-3, 0.034993],
                                   [3.53e-4, 6.872e-3, 0.028256]])

    extrapolation = copy.deepcopy(II)
    extrapolation_uncertainty = copy.deepcopy(dI)

    extrapolation[0,  0] = interp[0, 2] + interp[0, 1]*19. + interp[0, 0]*19.**2
    extrapolation[19, 0] = interp[0, 2]
    extrapolation[0,  1] = interp[1, 2] + interp[1, 1]*19. + interp[1, 0]*19.**2
    extrapolation[19, 1] = interp[1, 2]
    extrapolation[0,  2] = interp[2, 2] + interp[2, 1]*19. + interp[2, 0]*19.**2
    extrapolation[19, 2] = interp[2, 2]
    extrapolation[0, 3] = interp[3, 2] + interp[3, 1]*19. + interp[3, 0]*19.**2
    extrapolation[8, 3] = interp[3, 2] + interp[3, 1]*11. + interp[3, 0]*11.**2
    extrapolation[9, 3] = interp[3, 2] + interp[3, 1]*10. + interp[3, 0]*10.**2
    extrapolation[19, 3] = interp[3, 2]
    extrapolation[0, 4] = interp[4, 2] + interp[4, 1]*19. + interp[4, 0]*19.**2
    extrapolation[19, 4] = interp[4, 2]
    extrapolation[0, 5] = interp[5, 2] + interp[5, 1]*19. + interp[5, 0]*19.**2
    extrapolation[19, 5] = interp[5, 2]
    extrapolation[0, 6] = interp[6, 2] + interp[6, 1]*19. + interp[6, 0]*19.**2
    extrapolation[19, 6] = interp[6, 2]
    extrapolation[0, 7] = interp[7, 2] + interp[7, 1]*19. + interp[7, 0]*19.**2
    extrapolation[19, 7] = interp[7, 2]

    extrapolation_uncertainty[0, 0] = np.sqrt(interp_uncertainty[0, 2]**2 + (interp_uncertainty[0, 1]*19.)**2 +
                                              (interp_uncertainty[0, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 0] = np.sqrt(interp_uncertainty[0, 2]**2)
    extrapolation_uncertainty[0, 1] = np.sqrt(interp_uncertainty[1, 2]**2 + (interp_uncertainty[1, 1]*19.)**2 +
                                              (interp_uncertainty[1, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 1] = np.sqrt(interp_uncertainty[1, 2]**2)
    extrapolation_uncertainty[0, 2] = np.sqrt(interp_uncertainty[2, 2]**2 + (interp_uncertainty[2, 1]*19.)**2 +
                                              (interp_uncertainty[2, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 2] = np.sqrt(interp_uncertainty[2, 2]**2)
    extrapolation_uncertainty[0, 3] = np.sqrt(interp_uncertainty[3, 2]**2 + (interp_uncertainty[3, 1]*19.)**2 +
                                              (interp_uncertainty[3, 0]*19.**2)**2)
    extrapolation_uncertainty[8, 3] = np.sqrt(interp_uncertainty[3, 2]**2 + (interp_uncertainty[3, 1]*11.)**2 +
                                              (interp_uncertainty[3, 0]*11.**2)**2)
    extrapolation_uncertainty[9, 3] = np.sqrt(interp_uncertainty[3, 2]**2 + (interp_uncertainty[3, 1]*10.)**2 +
                                              (interp_uncertainty[3, 0]*10.**2)**2)
    extrapolation_uncertainty[19, 3] = np.sqrt(interp_uncertainty[3, 2]**2)
    extrapolation_uncertainty[0, 4] = np.sqrt(interp_uncertainty[4, 2]**2 + (interp_uncertainty[4, 1]*19.)**2 +
                                              (interp_uncertainty[4, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 4] = np.sqrt(interp_uncertainty[4, 2]**2)
    extrapolation_uncertainty[0, 5] = np.sqrt(interp_uncertainty[5, 2]**2 + (interp_uncertainty[5, 1]*19.)**2 +
                                              (interp_uncertainty[5, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 5] = np.sqrt(interp_uncertainty[5, 2]**2)
    extrapolation_uncertainty[0, 6] = np.sqrt(interp_uncertainty[6, 2]**2 + (interp_uncertainty[6, 1]*19.)**2 +
                                              (interp_uncertainty[6, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 6] = np.sqrt(interp_uncertainty[6, 2]**2)
    extrapolation_uncertainty[0, 7] = np.sqrt(interp_uncertainty[7, 2]**2 + (interp_uncertainty[7, 1]*19.)**2 +
                                              (interp_uncertainty[7, 0]*19.**2)**2)
    extrapolation_uncertainty[19, 7] = np.sqrt(interp_uncertainty[7, 2]**2)

    extrapolation[6, 1] = np.nan
    extrapolation[12, 5] = np.nan
    extrapolation_uncertainty[6, 1] = np.nan
    extrapolation_uncertainty[12, 5] = np.nan

    final_sensitivity = np.nanmean(extrapolation)
    n_elements = np.sum(np.logical_not(np.isnan(extrapolation_uncertainty)))
    final_sensitivity_uncertainty = np.sqrt(np.nansum(np.power(extrapolation_uncertainty, 2)))/n_elements
    result = extrapolation/final_sensitivity
    result_uncertainty = result * np.sqrt(np.square(extrapolation_uncertainty/extrapolation) +
                                          np.square(final_sensitivity_uncertainty/final_sensitivity))

    ws = workspace_with_instrument(axis_values=[1.,2.], intensities=ffm_with_mask,
                              uncertainties=ffm_uncertainty_with_mask, view='array')
    out = prepare_sensitivity_correction(ws, min_threshold=0.5, max_threshold=2.0)

    out_result = np.flip(np.transpose(out.extractY().reshape(8, 20)), 0)
    out_uncertainty = np.flip(np.transpose(out.extractE().reshape(8, 20)), 0)

    assert_allclose(result, out_result, equal_nan=True, atol=0.001)
    assert_allclose(result_uncertainty, out_uncertainty, equal_nan=True, atol=0.001)
