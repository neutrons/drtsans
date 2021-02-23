# Test drtsans.tof.eqsans.incoherence_correction_2d
import pytest
from drtsans.dataobjects import IQazimuthal
import numpy as np

import drtsans.tof.eqsans.incoherence_correction_2d as ic2d


def generate_test_data():
    # Generate test data
    # Created from
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/375b630f348d152cd18d919587bafd2a/test_inelastic_incoherent_avg_noerrestimation.xlsx
    # using np_array[qx, qy, lambda].flatten()
    # original data can be recreated by performing
    # i_of_q.intensity.reshape((11, 11, 5))
    intensity_vec = np.array([
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, 0.14, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, 0.13, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
        0.1, np.nan, np.nan, np.nan, np.nan,
    ])

    error_vec = np.sqrt(intensity_vec)

    vec_qx = np.arange(-0.1, 0.12, 0.02).repeat(55)
    vec_qy = np.tile(np.arange(-0.1, 0.12, 0.02).repeat(5), 11)

    wavelength_vec = np.tile(np.arange(1, 6, 1), 121)

    i_of_q = IQazimuthal(
        intensity=intensity_vec,
        error=error_vec,
        qx=vec_qx,
        qy=vec_qy,
        wavelength=wavelength_vec
    )

    return i_of_q


def test_gen_q_subset_mask():
    # q_subset of 2d case is each qx, qy where all lambda(qx, qy) is finite
    # in original test data, this can be represented as the ndarray ranges
    # 2:-2, 2, :
    # 2:-2, 8, :
    # 2, 2:-2, :
    # 8, 2:-2, :
    i_of_q = generate_test_data()
    q_subset_mask = ic2d.gen_q_subset_mask(i_of_q)
    # return mask to test data arangment for testing
    test_q_subset_mask = q_subset_mask.reshape((11, 11, 5))
    # Positive assertions
    assert test_q_subset_mask[2:-2, 2, :].all()
    assert test_q_subset_mask[2:-2, 8, :].all()
    assert test_q_subset_mask[2, 2:-2, :].all()
    assert test_q_subset_mask[8, 2:-2, :].all()
    # Negative assertions
    assert not test_q_subset_mask[:2, :, :].any()
    assert not test_q_subset_mask[-2:, :, :].any()
    assert not test_q_subset_mask[:, :2, :].any()
    assert not test_q_subset_mask[:, -2:, :].any()
    assert not test_q_subset_mask[3:8, 3:8, :].any()
    # Test filtering
    assert np.all(np.isfinite(i_of_q.intensity[q_subset_mask]))
    assert np.any(np.isnan(i_of_q.intensity[~q_subset_mask]))
