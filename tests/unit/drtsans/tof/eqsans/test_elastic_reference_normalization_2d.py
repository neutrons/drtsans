import pytest
from drtsans.dataobjects import IQazimuthal
from drtsans.tof.eqsans.elastic_reference_normalization import (
    normalize_by_elastic_reference2D,
)
import numpy as np


def create_testing_iq2d():
    """Create a test data I(Q, wavelength) as the attached EXCEL spreadsheet attached in gitlab story SANS 834
    https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/d268b5ddc440becf9c677e5e0e69e9b8/test_elastic_1_.xlsx
    Returns
    -------

    """
    # Elastic reference normalization vector
    k_vec = [1.0, 0.75, 1.4]
    k_error_vec = [0.0, 2.0e-4, 3.0e-4]

    # Intensity vector
    intensity_vec_3 = np.ones((11, 11), dtype=np.float64)
    intensity_vec_3[3:8, 3:8] = np.nan
    intensity_vec_4 = np.full((11, 11), np.nan, dtype=np.float64)
    intensity_vec_4[1:10, 1:10] = 1.2
    intensity_vec_4[2, 1:10] = 1.8
    intensity_vec_4[4:7, 4:7] = np.nan
    intensity_vec_5 = np.full((11, 11), np.nan, dtype=np.float64)
    intensity_vec_5[2:9, 2:9] = 0.9
    intensity_vec_5[5, 5] = np.nan
    intensity_vec = intensity_vec_3
    intensity_vec = np.concatenate((intensity_vec, intensity_vec_4), axis=1)
    intensity_vec = np.concatenate((intensity_vec, intensity_vec_5), axis=1)

    # Expected normalized intensity vector
    normalized_intensity_vec_4 = k_vec[1] * intensity_vec_4
    normalized_intensity_vec_5 = k_vec[2] * intensity_vec_5
    normalized_intensity_vec = intensity_vec_3
    normalized_intensity_vec = np.concatenate((normalized_intensity_vec, normalized_intensity_vec_4), axis=1)
    normalized_intensity_vec = np.concatenate((normalized_intensity_vec, normalized_intensity_vec_5), axis=1)

    # Error
    error_vec_3 = np.full((11, 11), 0.02, dtype=np.float64)
    error_vec_3[3:8, 3:8] = np.nan
    error_vec_4 = np.full((11, 11), np.nan, dtype=np.float64)
    error_vec_4[1:10, 1:10] = 0.03
    error_vec_4[4:7, 4:7] = np.nan
    error_vec_5 = np.full((11, 11), np.nan, dtype=np.float64)
    error_vec_5[2:9, 2:9] = 0.02
    error_vec_5[5, 5] = np.nan

    error_vec = error_vec_3
    error_vec = np.concatenate((error_vec, error_vec_4), axis=1)
    error_vec = np.concatenate((error_vec, error_vec_5), axis=1)

    # Expected normalized intensity error
    normalized_error_vec_3 = error_vec_3
    normalized_error_vec_4 = np.full((11, 11), np.nan, dtype=np.float64)
    normalized_error_vec_4[1:10, 1:10] = 0.0225015
    normalized_error_vec_4[2, 1:10] = 0.0225029
    normalized_error_vec_4[4:7, 4:7] = np.nan
    normalized_error_vec_5 = np.full((11, 11), np.nan, dtype=np.float64)
    normalized_error_vec_5[2:9, 2:9] = 0.0280013
    normalized_error_vec_5[5, 5] = np.nan
    normalized_error_vec = normalized_error_vec_3
    normalized_error_vec = np.concatenate((normalized_error_vec, normalized_error_vec_4), axis=1)
    normalized_error_vec = np.concatenate((normalized_error_vec, normalized_error_vec_5), axis=1)

    # Q vector
    vec_q = np.linspace(-0.1, 0.1, num=11, endpoint=True, dtype=np.float64)
    qx_matrix, qy_matrix = np.meshgrid(vec_q, vec_q, indexing="ij")
    qx_vec = qx_matrix
    qx_vec = np.concatenate((qx_vec, qx_matrix), axis=1)
    qx_vec = np.concatenate((qx_vec, qx_matrix), axis=1)
    qy_vec = qx_matrix
    qy_vec = np.concatenate((qy_vec, qy_matrix), axis=1)
    qy_vec = np.concatenate((qy_vec, qy_matrix), axis=1)

    # Wavelength vector
    wavelength_vec = np.full((11, 11), 3.0, dtype=np.float64)
    wavelength_vec = np.concatenate((wavelength_vec, np.full((11, 11), 4.0, dtype=np.float64)), axis=1)
    wavelength_vec = np.concatenate((wavelength_vec, np.full((11, 11), 5.0, dtype=np.float64)), axis=1)

    # Construct IQazimuthal
    i_of_q = IQazimuthal(intensity=intensity_vec, error=error_vec, qx=qx_vec, qy=qy_vec, wavelength=wavelength_vec)

    # Construct normalized IQazimuthal
    normalized_i_of_q = IQazimuthal(
        intensity=normalized_intensity_vec, error=normalized_error_vec, qx=qx_vec, qy=qy_vec, wavelength=wavelength_vec
    )

    return i_of_q, normalized_i_of_q, k_vec, k_error_vec


def test_normalize_by_elastic_reference2d():
    # Test normalizing I(Q2D) by elastic reference normalization factor K
    (
        i_of_q,
        expected_normalized_i_of_q,
        k_vec,
        k_error_vec,
    ) = create_testing_iq2d()  # noqa E501

    normalized_i_of_q = normalize_by_elastic_reference2D(i_of_q, k_vec, k_error_vec)

    np.testing.assert_allclose(normalized_i_of_q.intensity, expected_normalized_i_of_q.intensity, rtol=1e-5)
    np.testing.assert_allclose(normalized_i_of_q.error, expected_normalized_i_of_q.error, rtol=1e-5)


if __name__ == "__main__":
    pytest.main([__file__])
