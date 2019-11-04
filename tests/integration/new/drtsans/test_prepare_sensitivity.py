import numpy as np
import copy
from mantid.simpleapi import MaskDetectors
import pytest
from drtsans.sensitivity import calculate_sensitivity_correction


@pytest.mark.parametrize('workspace_with_instrument',
                          [dict(name='EQSANS', Nx=8, Ny=8)], indirect=True)
def test_prepare_sensitivity(workspace_with_instrument):
    flood_field_measurement = np.array([[5, 6, 7, 4, 6, 6,  5, 3],
                                        [7, 4, 6, 5, 5, 5,  4, 5],
                                        [3, 5, 7, 4, 6, 6,  5, 4],
                                        [6, 7, 5, 5, 4, 7,  4, 5],
                                        [3, 5, 7, 6, 7, 5, 56, 6],
                                        [7, 6, 5, 5, 6, 6,  4, 5],
                                        [4, 5, 4, 6, 5, 5,  5, 7],
                                        [5, 7, 6, 5, 6, 1,  4, 6]])
    flood_field_measurement_uncertainty = np.sqrt(flood_field_measurement)

    mask = np.ones((8, 8))
    mask[0, :] = np.nan
    mask[3, 3] = np.nan
    mask[4, 3] = np.nan
    mask[7, :] = np.nan

    ffm_with_mask = mask*flood_field_measurement
    ffm_uncertainty_with_mask = mask*flood_field_measurement_uncertainty
    F = np.nanmean(ffm_with_mask)
    n_elements = np.sum(np.logical_not(np.isnan(ffm_with_mask)))
    dF = np.sqrt(np.nansum(np.power(ffm_uncertainty_with_mask, 2)))/n_elements
    II = ffm_with_mask/F
    dI = II * np.sqrt(np.square(ffm_uncertainty_with_mask/ffm_with_mask) + np.square(dF/F))

    interp = np.array([[1.15043672, -0.26260662,  0.03137347],
                       [0.34048580,  0.28979423, -0.03532715],
                       [0.18591322,  0.27438438, -0.02357601],
                       [1.59838943, -0.39991784,  0.04023418],
                       [0.77215433,  0.05200745, -0.00730179],
                       [0.41601759,  0.23002406, -0.02466000],
                       [0.87065504, -0.08449276,  0.00819643],
                       [1.46153994, -0.24633847,  0.02046666]])

    interp_uncertainty = np.array([[1.01238129, 0.50205366, 0.05655452],
                                   [1.09729141, 0.53499255, 0.05800269],
                                   [1.04925415, 0.53275939, 0.05971851],
                                   [1.50183618, 0.79881347, 0.08807186],
                                   [1.08092149, 0.52927402, 0.05830330],
                                   [1.11487596, 0.54760600, 0.05991147],
                                   [1.10841878, 0.55543075, 0.06096038],
                                   [1.20327041, 0.57251127, 0.06194882]])

    extrapolation = copy.deepcopy(II)
    extrapolation_uncertainty = copy.deepcopy(dI)

    extrapolation[0, 0] = interp[0, 0] + interp[0, 1]*8. + interp[0, 2]*8.**2
    extrapolation[7, 0] = interp[0, 0] + interp[0, 1]*1. + interp[0, 2]*1.**2
    extrapolation[0, 1] = interp[1, 0] + interp[1, 1]*8. + interp[1, 2]*8.**2
    extrapolation[7, 1] = interp[1, 0] + interp[1, 1]*1. + interp[1, 2]*1.**2
    extrapolation[0, 2] = interp[2, 0] + interp[2, 1]*8. + interp[2, 2]*8.**2
    extrapolation[7, 2] = interp[2, 0] + interp[2, 1]*1. + interp[2, 2]*1.**2
    extrapolation[0, 3] = interp[3, 0] + interp[3, 1]*8. + interp[3, 2]*8.**2
    extrapolation[3, 3] = interp[3, 0] + interp[3, 1]*3. + interp[3, 2]*3.**2
    extrapolation[4, 3] = interp[3, 0] + interp[3, 1]*4. + interp[3, 2]*4.**2
    extrapolation[7, 3] = interp[3, 0] + interp[3, 1]*1. + interp[3, 2]*1.**2
    extrapolation[0, 4] = interp[4, 0] + interp[4, 1]*8. + interp[4, 2]*8.**2
    extrapolation[7, 4] = interp[4, 0] + interp[4, 1]*1. + interp[4, 2]*1.**2
    extrapolation[0, 5] = interp[5, 0] + interp[5, 1]*8. + interp[5, 2]*8.**2
    extrapolation[7, 5] = interp[5, 0] + interp[5, 1]*1. + interp[5, 2]*1.**2
    extrapolation[0, 6] = interp[6, 0] + interp[6, 1]*8. + interp[6, 2]*8.**2
    extrapolation[7, 6] = interp[6, 0] + interp[6, 1]*1. + interp[6, 2]*1.**2
    extrapolation[0, 7] = interp[7, 0] + interp[7, 1]*8. + interp[7, 2]*8.**2
    extrapolation[7, 7] = interp[7, 0] + interp[7, 1]*1. + interp[7, 2]*1.**2

    extrapolation[4, 6] = np.nan

    extrapolation_uncertainty[0, 0] = np.sqrt(interp_uncertainty[0, 0]**2 + (interp_uncertainty[0, 1]*8.)**2 +
                                              (interp_uncertainty[0, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 0] = np.sqrt(interp_uncertainty[0, 0]**2 + (interp_uncertainty[0, 1]*1.)**2 +
                                              (interp_uncertainty[0, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 1] = np.sqrt(interp_uncertainty[1, 0]**2 + (interp_uncertainty[1, 1]*8.)**2 +
                                              (interp_uncertainty[1, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 1] = np.sqrt(interp_uncertainty[1, 0]**2 + (interp_uncertainty[1, 1]*1.)**2 +
                                              (interp_uncertainty[1, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 2] = np.sqrt(interp_uncertainty[2, 0]**2 + (interp_uncertainty[2, 1]*8.)**2 +
                                              (interp_uncertainty[2, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 2] = np.sqrt(interp_uncertainty[2, 0]**2 + (interp_uncertainty[2, 1]*1.)**2 +
                                              (interp_uncertainty[2, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 3] = np.sqrt(interp_uncertainty[3, 0]**2 + (interp_uncertainty[3, 1]*8.)**2 +
                                              (interp_uncertainty[3, 2]*8.**2)**2)
    extrapolation_uncertainty[3, 3] = np.sqrt(interp_uncertainty[3, 0]**2 + (interp_uncertainty[3, 1]*3.)**2 +
                                              (interp_uncertainty[3, 2]*3.**2)**2)
    extrapolation_uncertainty[4, 3] = np.sqrt(interp_uncertainty[3, 0]**2 + (interp_uncertainty[3, 1]*4.)**2 +
                                              (interp_uncertainty[3, 2]*4.**2)**2)
    extrapolation_uncertainty[7, 3] = np.sqrt(interp_uncertainty[3, 0]**2 + (interp_uncertainty[3, 1]*1.)**2 +
                                              (interp_uncertainty[3, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 4] = np.sqrt(interp_uncertainty[4, 0]**2 + (interp_uncertainty[4, 1]*8.)**2 +
                                              (interp_uncertainty[4, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 4] = np.sqrt(interp_uncertainty[4, 0]**2 + (interp_uncertainty[4, 1]*1.)**2 +
                                              (interp_uncertainty[4, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 5] = np.sqrt(interp_uncertainty[5, 0]**2 + (interp_uncertainty[5, 1]*8.)**2 +
                                              (interp_uncertainty[5, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 5] = np.sqrt(interp_uncertainty[5, 0]**2 + (interp_uncertainty[5, 1]*1.)**2 +
                                              (interp_uncertainty[5, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 6] = np.sqrt(interp_uncertainty[6, 0]**2 + (interp_uncertainty[6, 1]*8.)**2 +
                                              (interp_uncertainty[6, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 6] = np.sqrt(interp_uncertainty[6, 0]**2 + (interp_uncertainty[6, 1]*1.)**2 +
                                              (interp_uncertainty[6, 2]*1.**2)**2)
    extrapolation_uncertainty[0, 7] = np.sqrt(interp_uncertainty[7, 0]**2 + (interp_uncertainty[7, 1]*8.)**2 +
                                              (interp_uncertainty[7, 2]*8.**2)**2)
    extrapolation_uncertainty[7, 7] = np.sqrt(interp_uncertainty[7, 0]**2 + (interp_uncertainty[7, 1]*1.)**2 +
                                              (interp_uncertainty[7, 2]*1.**2)**2)

    extrapolation_uncertainty[4, 6] = np.nan

    final_sensitivity = np.nanmean(extrapolation)
    n_elements = np.sum(np.logical_not(np.isnan(extrapolation_uncertainty)))
    final_sensitivity_uncertainty = np.sqrt(np.nansum(np.power(extrapolation_uncertainty, 2)))/n_elements
    result = extrapolation/final_sensitivity
    result_uncertainty = result * np.sqrt(np.square(extrapolation_uncertainty/extrapolation) +
                                          np.square(final_sensitivity_uncertainty/final_sensitivity))

    ws = workspace_with_instrument(axis_values=[1.,2.], intensities=ffm_with_mask,
                              uncertainties=ffm_uncertainty_with_mask, view='pixel')

    y = ws.extractY().flatten()
    #indices_to_mask = []
    #for i, yi in enumerate(y):
    #   if np.isnan(yi):
    #        indices_to_mask.append(i)
    indices_to_mask = np.arange(len(y))[np.isnan(y)]
    #print(y)
    MaskDetectors(ws, WorkspaceIndexList=indices_to_mask)
    np.set_printoptions(precision=4)
    #print(ws.extractY().reshape(8, 8))
    out = calculate_sensitivity_correction(ws, min_threashold=0.5, max_threshold=1.5)
    print(out.extractY().reshape(8, 8))
    print(out.extractE().reshape(8, 8))
    print(result)
    print(result_uncertainty)
    assert False
