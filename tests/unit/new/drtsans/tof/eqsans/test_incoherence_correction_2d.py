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
    intensity_vec = np.array(
        [
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            np.nan,
            np.nan,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            np.nan,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            0.14,
            0.11,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            0.14,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            0.13,
            0.15,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            0.13,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    error_vec = np.sqrt(intensity_vec)

    vec_qx = np.arange(-0.1, 0.12, 0.02).repeat(55)
    vec_qy = np.tile(np.arange(-0.1, 0.12, 0.02).repeat(5), 11)

    wavelength_vec = np.tile(np.arange(1, 6, 1), 121)

    i_of_q = IQazimuthal(
        intensity=intensity_vec,
        error=error_vec,
        qx=vec_qx,
        qy=vec_qy,
        wavelength=wavelength_vec,
    )

    return i_of_q


def test_reshape_q_azimuthal():
    # IQazimuthal does not require specific ordering; however,
    # calculations performed assume a specific ordering
    i_of_q = generate_test_data()
    # create some random order that will be de-randomed
    r_order = np.arange(i_of_q.intensity.shape[0])
    np.random.shuffle(r_order)
    r_i_of_q = IQazimuthal(
        intensity=i_of_q.intensity[r_order],
        error=i_of_q.error[r_order],
        qx=i_of_q.qx[r_order],
        qy=i_of_q.qy[r_order],
        wavelength=i_of_q.wavelength[r_order],
    )
    test_i_of_q = ic2d.reshape_q_azimuthal(r_i_of_q)
    i_drop_nan = i_of_q.intensity[np.isfinite(i_of_q.intensity)]
    test_i_drop_nan = test_i_of_q.intensity[np.isfinite(test_i_of_q.intensity)]
    assert np.array_equal(i_drop_nan, test_i_drop_nan)
    assert np.array_equal(i_of_q.qx, test_i_of_q.qx)
    assert np.array_equal(i_of_q.qy, test_i_of_q.qy)
    assert np.array_equal(i_of_q.wavelength, test_i_of_q.wavelength)


def test_gen_q_subset_mask():
    # q_subset of 2d case is each qx, qy where all lambda(qx, qy) is finite
    # in original test data, this can be represented as the ndarray ranges
    # 2:-2, 2, :
    # 2:-2, 8, :
    # 2, 2:-2, :
    # 8, 2:-2, :
    i_of_q = generate_test_data()
    qx_len = np.unique(i_of_q.qx).shape[0]
    qy_len = np.unique(i_of_q.qy).shape[0]
    wavelength_len = np.unique(i_of_q.wavelength).shape[0]
    q_subset_mask = ic2d.gen_q_subset_mask(i_of_q, qx_len, qy_len, wavelength_len)
    # return mask to test data arrangement for testing
    test_q_subset_mask = q_subset_mask.repeat(wavelength_len).reshape((11, 11, 5))
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
    q_subset_filter = q_subset_mask.repeat(wavelength_len)
    assert np.all(np.isfinite(i_of_q.intensity[q_subset_filter]))
    assert np.any(np.isnan(i_of_q.intensity[~q_subset_filter]))


def test_calculate_b2d():
    i_of_q = generate_test_data()
    qx_len = np.unique(i_of_q.qx).shape[0]
    qy_len = np.unique(i_of_q.qy).shape[0]
    wavelength_len = np.unique(i_of_q.wavelength).shape[0]
    q_subset_mask = ic2d.gen_q_subset_mask(i_of_q, qx_len, qy_len, wavelength_len)
    b_pack = ic2d.calculate_b2d(
        i_of_q, q_subset_mask, qx_len, qy_len, wavelength_len, min_incoh=False
    )
    b_vals, b_e_vals, ref = b_pack
    assert b_vals.shape[0] == wavelength_len
    assert b_e_vals.shape[0] == wavelength_len
    assert ref == 0
    known_b_vals = np.array([0, 0.03, 0.05, 0.04, 0.01])
    assert np.allclose(b_vals, known_b_vals)
    b_pack = ic2d.calculate_b2d(
        i_of_q, q_subset_mask, qx_len, qy_len, wavelength_len, min_incoh=True
    )
    b_vals, b_e_vals, ref = b_pack
    assert b_vals.shape[0] == wavelength_len
    assert b_e_vals.shape[0] == wavelength_len
    assert np.allclose(b_vals, known_b_vals)


if __name__ == "__main__":
    pytest.main([__file__])
