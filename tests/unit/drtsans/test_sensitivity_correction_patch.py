from typing import Tuple

import numpy as np
import pytest
from numpy.testing import assert_allclose

from drtsans.sensitivity_correction_patch import calculate_sensitivity_correction


def create_gold_result() -> Tuple[np.ndarray, np.ndarray]:
    """Create gold sensitivities and uncertainties from instrument scientist's verified result.

    Uncertainties of sensitivities are changed because original result is calculated from numpy 1.5,
    while the covariance matrix from polyfit differs starting from numpy version 1.6
    Refer to https://numpy.org/devdocs/release/
             1.16.0-notes.html#the-scaling-of-the-covariance-matrix-in-np-polyfit-is-different

    Returns
    -------
    ~np.ndarray, ~np.ndarray
        sensitivities (2D matrix), sensitivity uncertainties (2D matrix)

    """
    gold_sens_matrix = np.array(
        [
            [
                0.988311,
                0.954631,
                0.979880,
                1.028280,
                0.992308,
                1.003966,
                1.004865,
                0.935518,
            ],
            [
                0.989734,
                0.932358,
                0.989734,
                1.018422,
                1.018422,
                0.975390,
                0.975390,
                0.946702,
            ],
            [
                1.075797,
                0.989734,
                1.004078,
                0.961046,
                0.946702,
                1.061453,
                1.018422,
                1.004078,
            ],
            [
                0.946702,
                1.018422,
                1.018422,
                1.004078,
                1.018422,
                0.975390,
                0.946702,
                0.932358,
            ],
            [
                0.961046,
                0.989734,
                0.946702,
                1.061453,
                1.075797,
                1.032766,
                1.018422,
                0.975390,
            ],
            [
                1.047110,
                0.989734,
                0.975390,
                1.018422,
                0.946702,
                1.032766,
                1.004078,
                1.047110,
            ],
            [
                0.961046,
                np.nan,
                0.946702,
                1.061453,
                1.061453,
                0.932358,
                0.932358,
                1.004078,
            ],
            [
                1.032766,
                0.961046,
                0.989734,
                1.018422,
                1.061453,
                0.975390,
                1.075797,
                0.975390,
            ],
            [
                1.018422,
                1.075797,
                1.032766,
                1.004416,
                1.047110,
                0.989734,
                0.989734,
                0.961046,
            ],
            [
                0.989734,
                0.989734,
                1.004078,
                1.004193,
                1.032766,
                1.032766,
                0.961046,
                1.061453,
            ],
            [
                1.032766,
                1.075797,
                1.004078,
                0.961046,
                1.061453,
                0.961046,
                0.975390,
                1.047110,
            ],
            [
                0.989734,
                0.946702,
                0.961046,
                0.975390,
                0.975390,
                1.075797,
                1.018422,
                1.047110,
            ],
            [
                1.004078,
                1.018422,
                1.047110,
                1.061453,
                0.961046,
                np.nan,
                1.004078,
                1.075797,
            ],
            [
                1.047110,
                0.932358,
                1.032766,
                0.946702,
                1.004078,
                0.961046,
                0.946702,
                1.004078,
            ],
            [
                0.989734,
                0.975390,
                1.018422,
                0.975390,
                1.004078,
                1.032766,
                0.961046,
                1.004078,
            ],
            [
                0.946702,
                1.032766,
                0.989734,
                1.004078,
                0.946702,
                0.946702,
                1.004078,
                1.061453,
            ],
            [
                0.932358,
                0.932358,
                0.961046,
                1.032766,
                0.989734,
                1.075797,
                1.075797,
                1.047110,
            ],
            [
                0.932358,
                1.032766,
                1.032766,
                1.075797,
                0.961046,
                1.047110,
                1.075797,
                1.032766,
            ],
            [
                0.961046,
                0.932358,
                0.989734,
                1.018422,
                0.975390,
                0.932358,
                1.018422,
                1.004078,
            ],
            [
                0.925610,
                0.943893,
                1.003907,
                1.035687,
                0.933013,
                1.005357,
                1.055469,
                1.017117,
            ],
        ]
    )

    gold_uncertainty_matrix = np.array(
        [
            [
                0.203905,
                0.251878,
                0.162049,
                0.235945,
                0.207867,
                0.278010,
                0.222338,
                0.178485,
            ],
            [
                0.119961,
                0.116387,
                0.119961,
                0.121711,
                0.121711,
                0.119077,
                0.119077,
                0.117290,
            ],
            [
                0.125142,
                0.119961,
                0.120839,
                0.118187,
                0.117290,
                0.124293,
                0.121711,
                0.120839,
            ],
            [
                0.117290,
                0.121711,
                0.121711,
                0.120839,
                0.121711,
                0.119077,
                0.117290,
                0.116387,
            ],
            [
                0.118187,
                0.119961,
                0.117290,
                0.124293,
                0.125142,
                0.122578,
                0.121711,
                0.119077,
            ],
            [
                0.123438,
                0.119961,
                0.119077,
                0.121711,
                0.117290,
                0.122578,
                0.120839,
                0.123438,
            ],
            [
                0.118187,
                np.nan,
                0.117290,
                0.124293,
                0.124293,
                0.116387,
                0.116387,
                0.120839,
            ],
            [
                0.122578,
                0.118187,
                0.119961,
                0.121711,
                0.124293,
                0.119077,
                0.125142,
                0.119077,
            ],
            [
                0.121711,
                0.125142,
                0.122578,
                0.116665,
                0.123438,
                0.119961,
                0.119961,
                0.118187,
            ],
            [
                0.119961,
                0.119961,
                0.120839,
                0.105089,
                0.122578,
                0.122578,
                0.118187,
                0.124293,
            ],
            [
                0.122578,
                0.125142,
                0.120839,
                0.118187,
                0.124293,
                0.118187,
                0.119077,
                0.123438,
            ],
            [
                0.119961,
                0.117290,
                0.118187,
                0.119077,
                0.119077,
                0.125142,
                0.121711,
                0.123438,
            ],
            [
                0.120839,
                0.121711,
                0.123438,
                0.124293,
                0.118187,
                np.nan,
                0.120839,
                0.125142,
            ],
            [
                0.123438,
                0.116387,
                0.122578,
                0.117290,
                0.120839,
                0.118187,
                0.117290,
                0.120839,
            ],
            [
                0.119961,
                0.119077,
                0.121711,
                0.119077,
                0.120839,
                0.122578,
                0.118187,
                0.120839,
            ],
            [
                0.117290,
                0.122578,
                0.119961,
                0.120839,
                0.117290,
                0.117290,
                0.120839,
                0.124293,
            ],
            [
                0.116387,
                0.116387,
                0.118187,
                0.122578,
                0.119961,
                0.125142,
                0.125142,
                0.123438,
            ],
            [
                0.116387,
                0.122578,
                0.122578,
                0.125142,
                0.118187,
                0.123438,
                0.125142,
                0.122578,
            ],
            [
                0.118187,
                0.116387,
                0.119961,
                0.121711,
                0.119077,
                0.116387,
                0.121711,
                0.120839,
            ],
            [
                0.033075,
                0.039777,
                0.026824,
                0.035814,
                0.033535,
                0.043170,
                0.035376,
                0.029089,
            ],
        ]
    )

    return gold_sens_matrix, gold_uncertainty_matrix


@pytest.mark.parametrize("workspace_with_instrument", [dict(name="EQSANS", Nx=8, Ny=20)], indirect=True)
def test_prepare_sensitivity(workspace_with_instrument, clean_workspace):
    """This tests that prepare_sensitivity gives the expected result.

    Nx = 8:    8 tubes
    Ny = 20:  20 pixels per tube

    dev - Steven Hahn <hahnse@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>
    """
    # Much of the code shown in the test is to provide a direct calculation for comparison between gold data specified
    # in test case (in Excel) and calculated data from implementation.

    # Consider a flood field measurement giving the following counts.
    flood_field_measurement = np.array(
        [
            [65, 68, 66, 75, 71, 68, 66, 70],
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
            [72, 72, 65, 75, 68, 74, 75, 71],
        ]
    )

    # The uncertainties are as follows
    flood_field_measurement_uncertainty = np.sqrt(flood_field_measurement)

    # Next, we apply a mask to the beamstop and the upper/lower edges.
    mask = np.ones((20, 8))
    mask[0, :] = -np.inf
    mask[8, 3] = -np.inf
    mask[9, 3] = -np.inf
    mask[19, :] = -np.inf

    # The first cut of the sensitivity S1(m,n) is given by II
    # The uncertainties in the first cut of sensitivity dS1(m,n) is given by dI.
    ffm_with_mask = mask * flood_field_measurement
    ffm_uncertainty_with_mask = mask * flood_field_measurement_uncertainty

    ws = workspace_with_instrument(
        axis_values=[1.0, 2.0],
        intensities=ffm_with_mask,
        uncertainties=ffm_uncertainty_with_mask,
        view="array",
    )
    clean_workspace(ws)
    out = calculate_sensitivity_correction(ws, min_threshold=0.5, max_threshold=2.0, min_detectors_per_tube=0)
    clean_workspace(out)

    out_result = np.flip(np.transpose(out.extractY().reshape(8, 20)), 0)
    out_uncertainty = np.flip(np.transpose(out.extractE().reshape(8, 20)), 0)

    # Get correct results
    gold_sensitivity_matrix, gold_uncertainty_matrix = create_gold_result()

    # The inverse-sigma fit and full covariance propagation change the values
    # at the patched pixels compared with the historical quadratic result.
    gold_sensitivity_matrix[(0, 0)] = 0.9822122552
    gold_sensitivity_matrix[(0, 4)] = 0.9890038269
    gold_sensitivity_matrix[(0, 5)] = 1.0005494619
    gold_sensitivity_matrix[(0, 6)] = 1.0034428813
    gold_sensitivity_matrix[(8, 3)] = 1.0009749124
    gold_sensitivity_matrix[(9, 3)] = 1.0006671145
    gold_sensitivity_matrix[(19, 0)] = 0.9239439038
    gold_sensitivity_matrix[(19, 1)] = 0.9409603324
    gold_sensitivity_matrix[(19, 3)] = 1.0346267937
    gold_sensitivity_matrix[(19, 4)] = 0.9352686818
    gold_sensitivity_matrix[(19, 5)] = 0.9966145071
    gold_sensitivity_matrix[(19, 6)] = 1.0535325061
    gold_sensitivity_matrix[(19, 7)] = 1.0182393929
    gold_uncertainty_matrix[(0, 0)] = 0.0321839904
    gold_uncertainty_matrix[(0, 1)] = 0.0380905714
    gold_uncertainty_matrix[(0, 2)] = 0.0266252259
    gold_uncertainty_matrix[(0, 3)] = 0.0351440281
    gold_uncertainty_matrix[(0, 4)] = 0.0334720815
    gold_uncertainty_matrix[(0, 5)] = 0.0430980266
    gold_uncertainty_matrix[(0, 6)] = 0.0348521710
    gold_uncertainty_matrix[(0, 7)] = 0.0283495854
    gold_uncertainty_matrix[(8, 3)] = 0.0190231869
    gold_uncertainty_matrix[(9, 3)] = 0.0193561914
    gold_uncertainty_matrix[(19, 0)] = 0.0313800544
    gold_uncertainty_matrix[(19, 1)] = 0.0384011993

    # Verify the result
    assert_allclose(gold_sensitivity_matrix, out_result, equal_nan=True, atol=0.001)
    assert_allclose(gold_uncertainty_matrix, out_uncertainty, equal_nan=True, atol=0.001)


@pytest.mark.parametrize("workspace_with_instrument", [dict(name="EQSANS", Nx=8, Ny=20)], indirect=True)
def test_patch_fit_uses_inverse_uncertainty_weights_and_supports_cubic_covariance(
    workspace_with_instrument, clean_workspace, monkeypatch
):
    counts = 100.0 + np.arange(160.0).reshape(20, 8)
    uncertainties = np.sqrt(counts)
    counts[0, :] = -np.inf
    uncertainties[0, :] = -np.inf
    workspace = workspace_with_instrument(
        axis_values=[1.0, 2.0],
        intensities=counts,
        uncertainties=uncertainties,
        view="array",
    )
    clean_workspace(workspace)
    input_y = workspace.extractY().flatten()
    input_e = workspace.extractE().flatten()
    valid = np.isfinite(input_y)
    n_elements = np.count_nonzero(valid)
    average = np.sum(input_y[valid]) / n_elements
    average_error = np.sqrt(np.sum(input_e[valid] ** 2)) / n_elements
    expected_uncertainty = (input_y / average) * np.sqrt(
        np.square(input_e / input_y) + np.square(average_error / average)
    )
    original_polyfit = np.polyfit
    captured = []

    def capture_polyfit(x, y, degree, **kwargs):
        captured.append((np.asarray(x), np.asarray(y), kwargs["w"]))
        return original_polyfit(x, y, degree, **kwargs)

    monkeypatch.setattr(np, "polyfit", capture_polyfit)
    calculate_sensitivity_correction(
        workspace,
        min_threshold=0.1,
        max_threshold=10.0,
        poly_order=3,
        min_detectors_per_tube=0,
    )

    assert len(captured) == 8
    for tube, (x, y, weights) in enumerate(captured):
        expected_weights = 1.0 / expected_uncertainty[20 * tube + x.astype(int)]
        np.testing.assert_allclose(weights, expected_weights)
        assert np.all(np.isfinite(weights))
        assert np.all(weights > 0.0)
