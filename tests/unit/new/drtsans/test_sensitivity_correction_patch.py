import numpy as np
import copy
import pytest
from drtsans.sensitivity_correction_patch import calculate_sensitivity_correction
from numpy.testing import assert_allclose


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='EQSANS', Nx=8, Ny=20)], indirect=True)
def test_prepare_sensitivity_prototype(workspace_with_instrument):
    """This tests that prepare_sensitivity gives the expected result.

    Nx = 8:    8 tubes
    Ny = 20:  20 pixels per tube

    dev - Steven Hahn <hahnse@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>
    """
    # Much of the code shown in the test is to provide a direct calculation for comparison between gold data specified
    # in test case (in Excel) and calculated data from implementation.

    # Consider a flood field measurement giving the following counts.
    flood_field_measurement = np.array([[65, 68, 66, 75, 71, 68, 66, 70],
                                        [69, 65, 69, 71, 71, 68, 68, 66],
                                        [75, 69, 70, 67, 66, 74, 71, 70],
                                        [66, 71, 71, 70, 71, 68, 66, 65],
                                        [67, 69, 66, 74, 75, 72, 71, 68],
                                        [73, 69, 68, 71, 66, 72, 70, 73],
                                        [67, 6, 66, 74, 74, 65, 65, 70],
                                        [72, 67, 69, 71, 74, 68, 75, 68],
                                        [71, 75, 72, 73, 73, 69, 69, 67],
                                        [69, 69, 70, 66, 72, 72, 67, 74],
                                        [72, 75, 70, 67, 74, 67, 68, 73],
                                        [69, 66, 67, 68, 68, 75, 71, 73],
                                        [70, 71, 73, 74, 67, 492, 70, 75],
                                        [73, 65, 72, 66, 70, 67, 66, 70],
                                        [69, 68, 71, 68, 70, 72, 67, 70],
                                        [66, 72, 69, 70, 66, 66, 70, 74],
                                        [65, 65, 67, 72, 69, 75, 75, 73],
                                        [65, 72, 72, 75, 67, 73, 75, 72],
                                        [67, 65, 69, 71, 68, 65, 71, 70],
                                        [72, 72, 65, 75, 68, 74, 75, 71]])

    # The uncertainties are as follows
    flood_field_measurement_uncertainty = np.sqrt(flood_field_measurement)

    # Next, we apply a mask to the beamstop and the upper/lower edges.
    mask = np.ones((20, 8))
    mask[0, :] = np.NINF
    mask[8, 3] = np.NINF
    mask[9, 3] = np.NINF
    mask[19, :] = np.NINF

    # The first cut of the sensitivity S1(m,n) is given by II
    # The uncertainties in the first cut of sensitivity dS1(m,n) is given by dI.
    ffm_with_mask = mask * flood_field_measurement
    ffm_uncertainty_with_mask = mask * flood_field_measurement_uncertainty

    n_elements = ffm_with_mask.shape[0] * ffm_with_mask.shape[1] \
                 - np.count_nonzero(np.isnan(ffm_with_mask)) - np.count_nonzero(np.isneginf(ffm_with_mask))
    F = np.sum(
        [value for value in ffm_with_mask.ravel() if not np.isnan(value) and not np.isneginf(value)]) / n_elements
    dF = np.sqrt(np.sum([value ** 2 for value in ffm_uncertainty_with_mask.ravel()
                         if not np.isnan(value) and not np.isneginf(value)])) / n_elements
    II = ffm_with_mask / F
    dI = II * np.sqrt(np.square(ffm_uncertainty_with_mask / ffm_with_mask) + np.square(dF / F))

    # Using numpy.polyfit() with a 2nd-degree polynomial, one finds the following coefficients and uncertainties.
    interp = np.array([[-5.55e-4, 1.3720e-2, 0.892143],
                       [-6.55e-4, 1.2996e-2, 0.909765],
                       [-8.9e-5, 4.72e-4, 0.967609],
                       [2.96e-4, -5.991e-3, 0.998240],
                       [-6.63e-4, 1.5604e-2, 0.899279],
                       [0.000000, -6.4e-5, 0.969006],
                       [4.34e-4, -1.0815e-2, 1.017307],
                       [-5.71e-4, 6.709e-3, 0.980341]])

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

    # We apply the thresholds to S1(m,n).  The masked pixels are set to NaN
    extrapolation[6, 1] = np.nan
    extrapolation[12, 5] = np.nan
    extrapolation_uncertainty[6, 1] = np.nan
    extrapolation_uncertainty[12, 5] = np.nan

    # The patch is applied to the results of the previous step to produce S2(m,n).
    extrapolation[0, 0] = interp[0, 2] + interp[0, 1] * 19. + interp[0, 0] * 19. ** 2
    extrapolation[19, 0] = interp[0, 2]
    extrapolation[0, 1] = interp[1, 2] + interp[1, 1] * 19. + interp[1, 0] * 19. ** 2
    extrapolation[19, 1] = interp[1, 2]
    extrapolation[0, 2] = interp[2, 2] + interp[2, 1] * 19. + interp[2, 0] * 19. ** 2
    extrapolation[19, 2] = interp[2, 2]
    extrapolation[0, 3] = interp[3, 2] + interp[3, 1] * 19. + interp[3, 0] * 19. ** 2
    extrapolation[8, 3] = interp[3, 2] + interp[3, 1] * 11. + interp[3, 0] * 11. ** 2
    extrapolation[9, 3] = interp[3, 2] + interp[3, 1] * 10. + interp[3, 0] * 10. ** 2
    extrapolation[19, 3] = interp[3, 2]
    extrapolation[0, 4] = interp[4, 2] + interp[4, 1] * 19. + interp[4, 0] * 19. ** 2
    extrapolation[19, 4] = interp[4, 2]
    extrapolation[0, 5] = interp[5, 2] + interp[5, 1] * 19. + interp[5, 0] * 19. ** 2
    extrapolation[19, 5] = interp[5, 2]
    extrapolation[0, 6] = interp[6, 2] + interp[6, 1] * 19. + interp[6, 0] * 19. ** 2
    extrapolation[19, 6] = interp[6, 2]
    extrapolation[0, 7] = interp[7, 2] + interp[7, 1] * 19. + interp[7, 0] * 19. ** 2
    extrapolation[19, 7] = interp[7, 2]

    # The associated uncertainties, dS2(m,n) are given by the following.
    extrapolation_uncertainty[0, 0] = np.sqrt(interp_uncertainty[0, 2] ** 2 + (interp_uncertainty[0, 1] * 19.) ** 2 +
                                              (interp_uncertainty[0, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 0] = np.sqrt(interp_uncertainty[0, 2] ** 2)
    extrapolation_uncertainty[0, 1] = np.sqrt(interp_uncertainty[1, 2] ** 2 + (interp_uncertainty[1, 1] * 19.) ** 2 +
                                              (interp_uncertainty[1, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 1] = np.sqrt(interp_uncertainty[1, 2] ** 2)
    extrapolation_uncertainty[0, 2] = np.sqrt(interp_uncertainty[2, 2] ** 2 + (interp_uncertainty[2, 1] * 19.) ** 2 +
                                              (interp_uncertainty[2, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 2] = np.sqrt(interp_uncertainty[2, 2] ** 2)
    extrapolation_uncertainty[0, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 19.) ** 2 +
                                              (interp_uncertainty[3, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[8, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 11.) ** 2 +
                                              (interp_uncertainty[3, 0] * 11. ** 2) ** 2)
    extrapolation_uncertainty[9, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 10.) ** 2 +
                                              (interp_uncertainty[3, 0] * 10. ** 2) ** 2)
    extrapolation_uncertainty[19, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2)
    extrapolation_uncertainty[0, 4] = np.sqrt(interp_uncertainty[4, 2] ** 2 + (interp_uncertainty[4, 1] * 19.) ** 2 +
                                              (interp_uncertainty[4, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 4] = np.sqrt(interp_uncertainty[4, 2] ** 2)
    extrapolation_uncertainty[0, 5] = np.sqrt(interp_uncertainty[5, 2] ** 2 + (interp_uncertainty[5, 1] * 19.) ** 2 +
                                              (interp_uncertainty[5, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 5] = np.sqrt(interp_uncertainty[5, 2] ** 2)
    extrapolation_uncertainty[0, 6] = np.sqrt(interp_uncertainty[6, 2] ** 2 + (interp_uncertainty[6, 1] * 19.) ** 2 +
                                              (interp_uncertainty[6, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 6] = np.sqrt(interp_uncertainty[6, 2] ** 2)
    extrapolation_uncertainty[0, 7] = np.sqrt(interp_uncertainty[7, 2] ** 2 + (interp_uncertainty[7, 1] * 19.) ** 2 +
                                              (interp_uncertainty[7, 0] * 19. ** 2) ** 2)
    extrapolation_uncertainty[19, 7] = np.sqrt(interp_uncertainty[7, 2] ** 2)

    # The final sensitivity, S(m,n), is produced by dividing this result
    # by the average value per Equations A3.13 and A3.14
    n_elements = ffm_with_mask.shape[0] * ffm_with_mask.shape[1] \
                 - np.count_nonzero(np.isnan(extrapolation)) - np.count_nonzero(np.isneginf(extrapolation))
    final_sensitivity = np.sum([value for value in extrapolation.ravel()
                                if not np.isnan(value) and not np.isneginf(value)]) / n_elements
    final_sensitivity_uncertainty = np.sqrt(np.sum([value ** 2 for value in extrapolation_uncertainty.ravel()
                                                    if not np.isnan(value) and not np.isneginf(value)])) / n_elements
    result = extrapolation / final_sensitivity
    result_uncertainty = result * np.sqrt(np.square(extrapolation_uncertainty / extrapolation) +
                                          np.square(final_sensitivity_uncertainty / final_sensitivity))
    ws = workspace_with_instrument(axis_values=[1., 2.], intensities=ffm_with_mask,
                                   uncertainties=ffm_uncertainty_with_mask, view='array')
    out = calculate_sensitivity_correction(ws, min_threshold=0.5, max_threshold=2.0,
                                           min_detectors_per_tube=0)

    out_result = np.flip(np.transpose(out.extractY().reshape(8, 20)), 0)
    out_uncertainty = np.flip(np.transpose(out.extractE().reshape(8, 20)), 0)

    assert_allclose(result, out_result, equal_nan=True, atol=0.001)
    assert_allclose(result_uncertainty, out_uncertainty, equal_nan=True, atol=0.001)


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='EQSANS', Nx=8, Ny=20)], indirect=True)
def test_prepare_sensitivity(workspace_with_instrument):
    """This tests that prepare_sensitivity gives the expected result.

    Nx = 8:    8 tubes
    Ny = 20:  20 pixels per tube

    dev - Steven Hahn <hahnse@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>
    """
    # Much of the code shown in the test is to provide a direct calculation for comparison between gold data specified
    # in test case (in Excel) and calculated data from implementation.

    # Consider a flood field measurement giving the following counts.
    flood_field_measurement = np.array([[65, 68, 66, 75, 71, 68, 66, 70],
                                        [69, 65, 69, 71, 71, 68, 68, 66],
                                        [75, 69, 70, 67, 66, 74, 71, 70],
                                        [66, 71, 71, 70, 71, 68, 66, 65],
                                        [67, 69, 66, 74, 75, 72, 71, 68],
                                        [73, 69, 68, 71, 66, 72, 70, 73],
                                        [67, 6, 66, 74, 74, 65, 65, 70],
                                        [72, 67, 69, 71, 74, 68, 75, 68],
                                        [71, 75, 72, 73, 73, 69, 69, 67],
                                        [69, 69, 70, 66, 72, 72, 67, 74],
                                        [72, 75, 70, 67, 74, 67, 68, 73],
                                        [69, 66, 67, 68, 68, 75, 71, 73],
                                        [70, 71, 73, 74, 67, 492, 70, 75],
                                        [73, 65, 72, 66, 70, 67, 66, 70],
                                        [69, 68, 71, 68, 70, 72, 67, 70],
                                        [66, 72, 69, 70, 66, 66, 70, 74],
                                        [65, 65, 67, 72, 69, 75, 75, 73],
                                        [65, 72, 72, 75, 67, 73, 75, 72],
                                        [67, 65, 69, 71, 68, 65, 71, 70],
                                        [72, 72, 65, 75, 68, 74, 75, 71]])

    # The uncertainties are as follows
    flood_field_measurement_uncertainty = np.sqrt(flood_field_measurement)

    # Next, we apply a mask to the beamstop and the upper/lower edges.
    mask = np.ones((20, 8))
    mask[0, :] = np.NINF
    mask[8, 3] = np.NINF
    mask[9, 3] = np.NINF
    mask[19, :] = np.NINF

    # The first cut of the sensitivity S1(m,n) is given by II
    # The uncertainties in the first cut of sensitivity dS1(m,n) is given by dI.
    ffm_with_mask = mask * flood_field_measurement
    ffm_uncertainty_with_mask = mask * flood_field_measurement_uncertainty

    # n_elements = ffm_with_mask.shape[0] * ffm_with_mask.shape[1] \
    #     - np.count_nonzero(np.isnan(ffm_with_mask)) - np.count_nonzero(np.isneginf(ffm_with_mask))
    # F = np.sum(
    #     [value for value in ffm_with_mask.ravel() if not np.isnan(value) and not np.isneginf(value)]) / n_elements
    # dF = np.sqrt(np.sum([value ** 2 for value in ffm_uncertainty_with_mask.ravel()
    #                      if not np.isnan(value) and not np.isneginf(value)])) / n_elements
    # II = ffm_with_mask / F
    # dI = II * np.sqrt(np.square(ffm_uncertainty_with_mask / ffm_with_mask) + np.square(dF / F))
    #
    # # Using numpy.polyfit() with a 2nd-degree polynomial, one finds the following coefficients and uncertainties.
    # interp = np.array([[-5.55e-4, 1.3720e-2, 0.892143],
    #                    [-6.55e-4, 1.2996e-2, 0.909765],
    #                    [-8.9e-5, 4.72e-4, 0.967609],
    #                    [2.96e-4, -5.991e-3, 0.998240],
    #                    [-6.63e-4, 1.5604e-2, 0.899279],
    #                    [0.000000, -6.4e-5, 0.969006],
    #                    [4.34e-4, -1.0815e-2, 1.017307],
    #                    [-5.71e-4, 6.709e-3, 0.980341]])
    #
    # interp_uncertainty = np.array([[4.01e-4, 7.882e-3, 0.032903],
    #                                [4.99e-4, 9.798e-3, 0.040252],
    #                                [3.19e-4, 6.248e-3, 0.025796],
    #                                [4.72e-4, 9.217e-3, 0.035956],
    #                                [4.10e-4, 8.021e-3, 0.033378],
    #                                [5.53e-4, 1.0765e-2, 0.043732],
    #                                [4.40e-4, 8.559e-3, 0.034993],
    #                                [3.53e-4, 6.872e-3, 0.028256]])
    #
    # extrapolation = copy.deepcopy(II)
    # extrapolation_uncertainty = copy.deepcopy(dI)
    #
    # # We apply the thresholds to S1(m,n).  The masked pixels are set to NaN
    # extrapolation[6, 1] = np.nan
    # extrapolation[12, 5] = np.nan
    # extrapolation_uncertainty[6, 1] = np.nan
    # extrapolation_uncertainty[12, 5] = np.nan
    #
    # # The patch is applied to the results of the previous step to produce S2(m,n).
    # extrapolation[0, 0] = interp[0, 2] + interp[0, 1] * 19. + interp[0, 0] * 19. ** 2
    # extrapolation[19, 0] = interp[0, 2]
    # extrapolation[0, 1] = interp[1, 2] + interp[1, 1] * 19. + interp[1, 0] * 19. ** 2
    # extrapolation[19, 1] = interp[1, 2]
    # extrapolation[0, 2] = interp[2, 2] + interp[2, 1] * 19. + interp[2, 0] * 19. ** 2
    # extrapolation[19, 2] = interp[2, 2]
    # extrapolation[0, 3] = interp[3, 2] + interp[3, 1] * 19. + interp[3, 0] * 19. ** 2
    # extrapolation[8, 3] = interp[3, 2] + interp[3, 1] * 11. + interp[3, 0] * 11. ** 2
    # extrapolation[9, 3] = interp[3, 2] + interp[3, 1] * 10. + interp[3, 0] * 10. ** 2
    # extrapolation[19, 3] = interp[3, 2]
    # extrapolation[0, 4] = interp[4, 2] + interp[4, 1] * 19. + interp[4, 0] * 19. ** 2
    # extrapolation[19, 4] = interp[4, 2]
    # extrapolation[0, 5] = interp[5, 2] + interp[5, 1] * 19. + interp[5, 0] * 19. ** 2
    # extrapolation[19, 5] = interp[5, 2]
    # extrapolation[0, 6] = interp[6, 2] + interp[6, 1] * 19. + interp[6, 0] * 19. ** 2
    # extrapolation[19, 6] = interp[6, 2]
    # extrapolation[0, 7] = interp[7, 2] + interp[7, 1] * 19. + interp[7, 0] * 19. ** 2
    # extrapolation[19, 7] = interp[7, 2]
    #
    # # The associated uncertainties, dS2(m,n) are given by the following.
    # extrapolation_uncertainty[0, 0] = np.sqrt(interp_uncertainty[0, 2] ** 2 + (interp_uncertainty[0, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[0, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 0] = np.sqrt(interp_uncertainty[0, 2] ** 2)
    # extrapolation_uncertainty[0, 1] = np.sqrt(interp_uncertainty[1, 2] ** 2 + (interp_uncertainty[1, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[1, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 1] = np.sqrt(interp_uncertainty[1, 2] ** 2)
    # extrapolation_uncertainty[0, 2] = np.sqrt(interp_uncertainty[2, 2] ** 2 + (interp_uncertainty[2, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[2, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 2] = np.sqrt(interp_uncertainty[2, 2] ** 2)
    # extrapolation_uncertainty[0, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[3, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[8, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 11.) ** 2 +
    #                                           (interp_uncertainty[3, 0] * 11. ** 2) ** 2)
    # extrapolation_uncertainty[9, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2 + (interp_uncertainty[3, 1] * 10.) ** 2 +
    #                                           (interp_uncertainty[3, 0] * 10. ** 2) ** 2)
    # extrapolation_uncertainty[19, 3] = np.sqrt(interp_uncertainty[3, 2] ** 2)
    # extrapolation_uncertainty[0, 4] = np.sqrt(interp_uncertainty[4, 2] ** 2 + (interp_uncertainty[4, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[4, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 4] = np.sqrt(interp_uncertainty[4, 2] ** 2)
    # extrapolation_uncertainty[0, 5] = np.sqrt(interp_uncertainty[5, 2] ** 2 + (interp_uncertainty[5, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[5, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 5] = np.sqrt(interp_uncertainty[5, 2] ** 2)
    # extrapolation_uncertainty[0, 6] = np.sqrt(interp_uncertainty[6, 2] ** 2 + (interp_uncertainty[6, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[6, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 6] = np.sqrt(interp_uncertainty[6, 2] ** 2)
    # extrapolation_uncertainty[0, 7] = np.sqrt(interp_uncertainty[7, 2] ** 2 + (interp_uncertainty[7, 1] * 19.) ** 2 +
    #                                           (interp_uncertainty[7, 0] * 19. ** 2) ** 2)
    # extrapolation_uncertainty[19, 7] = np.sqrt(interp_uncertainty[7, 2] ** 2)
    #
    # # The final sensitivity, S(m,n), is produced by dividing this result
    # # by the average value per Equations A3.13 and A3.14
    # n_elements = ffm_with_mask.shape[0] * ffm_with_mask.shape[1] \
    #     - np.count_nonzero(np.isnan(extrapolation)) - np.count_nonzero(np.isneginf(extrapolation))
    # final_sensitivity = np.sum([value for value in extrapolation.ravel()
    #                             if not np.isnan(value) and not np.isneginf(value)]) / n_elements
    # final_sensitivity_uncertainty = np.sqrt(np.sum([value ** 2 for value in extrapolation_uncertainty.ravel()
    #                                                 if not np.isnan(value) and not np.isneginf(value)])) / n_elements
    # result = extrapolation / final_sensitivity
    # result_uncertainty = result * np.sqrt(np.square(extrapolation_uncertainty / extrapolation) +
    #                                       np.square(final_sensitivity_uncertainty / final_sensitivity))
    ws = workspace_with_instrument(axis_values=[1., 2.], intensities=ffm_with_mask,
                                   uncertainties=ffm_uncertainty_with_mask, view='array')
    out = calculate_sensitivity_correction(ws, min_threshold=0.5, max_threshold=2.0,
                                           min_detectors_per_tube=0)

    out_result = np.flip(np.transpose(out.extractY().reshape(8, 20)), 0)
    out_uncertainty = np.flip(np.transpose(out.extractE().reshape(8, 20)), 0)

    print('Shape[Out] = {}, {}'.format(out_result.shape, out_uncertainty.shape()))

    # assert_allclose(result, out_result, equal_nan=True, atol=0.001)
    # assert_allclose(result_uncertainty, out_uncertainty, equal_nan=True, atol=0.001)

    assert 1 == 2
