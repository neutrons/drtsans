import os
import pytest
from drtsans.dataobjects import IQmod

from drtsans.tof.eqsans.elastic_reference_normalization import (
    determine_common_mod_q_range_mesh,
    normalize_by_elastic_reference,
)
from drtsans.tof.eqsans.elastic_reference_normalization import (
    calculate_scale_factor_mesh_grid,
)
from drtsans.tof.eqsans.elastic_reference_normalization import (
    reshape_q_wavelength_matrix,
)
from drtsans.tof.eqsans.elastic_reference_normalization import (
    determine_reference_wavelength_q1d_mesh,
)
from drtsans.tof.eqsans.elastic_reference_normalization import normalize_intensity_q1d
import numpy as np


def test_reshaped_imodq_2d():
    """Test the method to reshape I(Q) to 2D  mesh grid:"""
    # Get test data
    test_data = create_testing_iq1d()
    i_of_q = test_data[0]

    # Reshape
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_of_q)

    # Verify shape
    assert wl_vec.shape == (4,)
    assert q_vec.shape == (20,)
    assert i_array.shape == (20, 4)
    assert error_array.shape == (20, 4)

    # Verify value
    # intensity from q = 0.07 to 0.12
    gold_intensity_array = np.array(
        [
            [28.727, 31.600, 34.473, 25.855],
            [18.262, 20.088, 21.914, 16.435],
            [12.076, 13.283, 14.491, 10.868],
            [8.264, 9.091, 9.917, 7.438],
            [5.827, 6.410, 6.993, 5.244],
        ]
    )
    for i_q in range(6, 11):
        np.testing.assert_allclose(i_array[i_q], gold_intensity_array[i_q - 6])
    assert i_array[19, 0] == 0.595

    # error for wl = 3
    wl3_errors = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            5.360,
            4.273,
            3.475,
            2.875,
            2.414,
            2.053,
            1.767,
            1.535,
            1.346,
            1.189,
            1.058,
            0.947,
            0.852,
            0.771,
        ]
    )
    np.testing.assert_allclose(error_array[:, 0], wl3_errors, equal_nan=True)
    assert error_array[0, 3] == 27.273
    assert error_array[16, 1] == 1.028


def test_determine_common_q_range():
    """Test method to determine q range common to all wavelength"""
    # Get test data
    test_data = create_testing_iq1d()
    i_of_q = test_data[0]

    # Reshape
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_of_q)

    # Call the mesh-grid algorithm
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # Verify by comparing results from 2 algorithms
    assert q_vec[qmin_index] == 0.07
    assert q_vec[qmax_index] == 0.11


def test_determine_reference_wavelength():
    """Test method to determine the reference wavelength's intensity and error
    for each momentum transfer Q1D
    """
    # Get test data
    test_data = create_testing_iq1d()
    i_of_q = test_data[0]

    # Reshape
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(i_of_q)

    # Calculate reference wavelength and the related intensities and errors
    ref_wl_intensities = determine_reference_wavelength_q1d_mesh(
        wl_vec, q_vec, i_array, error_array, qmin_index=6, qmax_index=10
    )

    # Verify
    # Must be all at wl = 3.0
    np.testing.assert_allclose(ref_wl_intensities.ref_wl_vec, np.zeros_like(q_vec) + 3.0)
    np.testing.assert_allclose(
        ref_wl_intensities.intensity_vec[6:11],
        np.array([28.727, 18.262, 12.076, 8.264, 5.827]),
    )


def test_calculate_scale_factor():
    """Test the method to calculate scale factor K(wavelength) and delta K(wavelength)"""
    # Get testing data and gold data
    test_i_of_q, gold_k_vec, gold_k_error_vec, gold_intensity_vec, gold_error_vec = create_testing_iq1d()
    # Reshape
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_i_of_q)

    # Calculate Qmin and Qmax
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # Calculate reference
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array, qmin_index, qmax_index)

    # Calculate scale factor
    k_vec, k_error_vec = calculate_scale_factor_mesh_grid(
        wl_vec, i_array, error_array, ref_wl_ie, qmin_index, qmax_index
    )

    np.testing.assert_allclose(k_vec, gold_k_vec, rtol=1e-5)
    np.testing.assert_allclose(k_error_vec, gold_k_error_vec, rtol=1e-5)


def test_workflow_q1d(temp_directory):
    """Test method normalize_by_elastic_reference"""
    # Get testing data and gold data
    test_i_of_q, gold_k_vec, gold_k_error_vec, gold_intensity_vec, gold_error_vec = create_testing_iq1d()

    output_dir = temp_directory()

    # Normalize
    normalized_iq1d = normalize_by_elastic_reference(
        test_i_of_q,
        gold_k_vec,
        gold_k_error_vec,
        output_wavelength_dependent_profile=True,
        output_dir=output_dir,
    )

    # Verify
    np.testing.assert_allclose(normalized_iq1d.mod_q, test_i_of_q.mod_q, rtol=1e-6)
    np.testing.assert_allclose(normalized_iq1d.wavelength, test_i_of_q.wavelength, rtol=1e-6)
    np.testing.assert_allclose(normalized_iq1d.intensity, gold_intensity_vec, 1e-3)
    np.testing.assert_allclose(normalized_iq1d.error, gold_error_vec, 1e-3)

    # check the output iq wavelength profiles
    expected_len = [14, 13, 13, 11]
    for n in range(4):
        wl = n + 3.0
        filename = os.path.join(output_dir, f"IQ_{wl:.3f}_before_k_correction.dat")
        assert os.path.exists(filename)
        data = np.loadtxt(filename)
        assert len(data) == expected_len[n]

        filename = os.path.join(output_dir, f"IQ_{wl:.3f}_after_k_correction.dat")
        assert os.path.exists(filename)
        data = np.loadtxt(filename)
        assert len(data) == expected_len[n]


def test_normalize_i_of_q1d():
    """Test the refined method to normalize I(Q1D)"""
    # Get testing data and gold data
    test_i_of_q, gold_k_vec, gold_k_error_vec, gold_intensity_vec, gold_error_vec = create_testing_iq1d()
    # Reshape
    wl_vec, q_vec, i_array, error_array, dq_array = reshape_q_wavelength_matrix(test_i_of_q)

    # Calculate Qmin and Qmax
    qmin_index, qmax_index = determine_common_mod_q_range_mesh(q_vec, i_array)

    # Calculate reference
    ref_wl_ie = determine_reference_wavelength_q1d_mesh(wl_vec, q_vec, i_array, error_array, qmin_index, qmax_index)

    # Calculate scale factor
    k_vec, k_error_vec = calculate_scale_factor_mesh_grid(
        wl_vec, i_array, error_array, ref_wl_ie, qmin_index, qmax_index
    )

    # Normalize
    normalized = normalize_intensity_q1d(
        wl_vec,
        q_vec,
        i_array,
        error_array,
        k_vec,
        k_error_vec,
    )
    normalized_intensity = normalized[0].flatten()
    normalized_error = normalized[1].flatten()

    np.testing.assert_allclose(normalized_intensity, gold_intensity_vec, rtol=8e-4, equal_nan=True)
    np.testing.assert_allclose(normalized_error, gold_error_vec, rtol=1e-3, equal_nan=True)


def create_testing_iq1d():
    """Create a test data I(Q, wavelength) as the attached EXCEL spreadsheet attached in gitlab story
    Returns
    -------

    """
    # Intensity vector
    intensity_vec = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            743.802,
            np.nan,
            np.nan,
            np.nan,
            459.184,
            np.nan,
            np.nan,
            332.410,
            249.307,
            np.nan,
            np.nan,
            177.515,
            133.136,
            np.nan,
            89.796,
            97.959,
            73.469,
            np.nan,
            51.985,
            56.711,
            42.533,
            28.727,
            31.600,
            34.473,
            25.855,
            18.262,
            20.088,
            21.914,
            16.435,
            12.076,
            13.283,
            14.491,
            10.868,
            8.264,
            9.091,
            9.917,
            7.438,
            5.827,
            6.410,
            6.993,
            5.244,
            4.217,
            4.638,
            5.060,
            np.nan,
            3.121,
            3.433,
            3.745,
            np.nan,
            2.356,
            2.592,
            2.828,
            np.nan,
            1.811,
            1.992,
            2.173,
            np.nan,
            1.413,
            1.555,
            np.nan,
            np.nan,
            1.119,
            1.230,
            np.nan,
            np.nan,
            0.896,
            np.nan,
            np.nan,
            np.nan,
            0.727,
            np.nan,
            np.nan,
            np.nan,
            0.595,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    print(f"intensity: {intensity_vec.shape}")

    # Error
    error_vec = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            27.273,
            np.nan,
            np.nan,
            np.nan,
            21.429,
            np.nan,
            np.nan,
            18.232,
            15.789,
            np.nan,
            np.nan,
            13.323,
            11.538,
            np.nan,
            16.792,
            9.897,
            8.571,
            np.nan,
            10.611,
            7.531,
            6.522,
            5.360,
            4.605,
            5.871,
            5.085,
            4.273,
            4.374,
            4.681,
            4.054,
            3.475,
            3.640,
            3.807,
            3.297,
            2.875,
            2.985,
            3.149,
            2.727,
            2.414,
            2.470,
            2.644,
            2.290,
            2.053,
            2.095,
            2.249,
            np.nan,
            1.767,
            1.772,
            1.935,
            np.nan,
            1.535,
            1.522,
            1.682,
            np.nan,
            1.346,
            1.322,
            1.474,
            np.nan,
            1.189,
            1.161,
            np.nan,
            np.nan,
            1.058,
            1.028,
            np.nan,
            np.nan,
            0.947,
            np.nan,
            np.nan,
            np.nan,
            0.852,
            np.nan,
            np.nan,
            np.nan,
            0.771,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    # Q vector
    vec_q = np.arange(1, 21) * 0.01
    vec_q = np.repeat(vec_q, 4)
    error_q_vec = np.sqrt(vec_q)
    print(f"q: {vec_q.shape}")

    # Wavelength vector
    wavelength_vec = np.arange(3, 7) * 1.0
    wavelength_vec = np.tile(wavelength_vec, 20)
    print(f"lambda: {wavelength_vec.shape}")

    # Construct IQmod
    i_of_q = IQmod(
        intensity=intensity_vec, error=error_vec, mod_q=vec_q, delta_mod_q=error_q_vec, wavelength=wavelength_vec
    )

    # Expected K vector
    expected_k_vec = np.array([1.0000, 0.9090909090909, 0.833333333333333, 1.11111111111111])

    expected_k_error_vec = np.array([0.18072656, 0.1506194, 0.14418833, 0.20631421])

    expected_corrected_intensities = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            826.446,
            np.nan,
            np.nan,
            np.nan,
            510.204,
            np.nan,
            np.nan,
            277.008,
            277.008,
            np.nan,
            np.nan,
            147.929,
            147.929,
            np.nan,
            81.633,
            81.633,
            81.633,
            np.nan,
            47.259,
            47.259,
            47.259,
            28.727,
            28.727,
            28.727,
            28.727,
            18.262,
            18.262,
            18.262,
            18.262,
            12.076,
            12.076,
            12.076,
            12.076,
            8.264,
            8.264,
            8.264,
            8.264,
            5.827,
            5.827,
            5.827,
            5.827,
            4.217,
            4.217,
            4.217,
            np.nan,
            3.121,
            3.121,
            3.121,
            np.nan,
            2.356,
            2.356,
            2.356,
            np.nan,
            1.811,
            1.811,
            1.811,
            np.nan,
            1.413,
            1.413,
            np.nan,
            np.nan,
            1.119,
            1.119,
            np.nan,
            np.nan,
            0.896,
            np.nan,
            np.nan,
            np.nan,
            0.727,
            np.nan,
            np.nan,
            np.nan,
            0.595,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    expected_corrected_errors = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            156.42031992045736,
            np.nan,
            np.nan,
            np.nan,
            97.68243958016596,
            np.nan,
            np.nan,
            50.280060338882606,
            54.34506514416749,
            np.nan,
            np.nan,
            27.89978362627843,
            30.312283927845428,
            np.nan,
            20.39507103108296,
            16.356129961569046,
            17.901100412261805,
            np.nan,
            12.424166260136452,
            10.307757852900577,
            11.38057418635211,
            5.36,
            6.338696661245461,
            6.974460773306193,
            7.770241519009199,
            4.273,
            4.996587808230468,
            5.019984875333805,
            5.638022465861475,
            3.475,
            3.866876044267073,
            3.798728164584181,
            4.295058472897778,
            2.875,
            3.0395247176703917,
            2.9884454307509265,
            3.39643309869897,
            2.414,
            2.44421158125286,
            2.423074887310439,
            2.7649064627527142,
            2.053,
            2.0286144409287212,
            2.011157337753124,
            np.nan,
            1.767,
            1.6918573372416295,
            1.70050032072628,
            np.nan,
            1.535,
            1.4376561036235684,
            1.4597643383392334,
            np.nan,
            1.346,
            1.238700522312652,
            1.2676554475094703,
            np.nan,
            1.189,
            1.0811261397078902,
            np.nan,
            np.nan,
            1.058,
            0.9527287998571962,
            np.nan,
            np.nan,
            0.947,
            np.nan,
            np.nan,
            np.nan,
            0.852,
            np.nan,
            np.nan,
            np.nan,
            0.771,
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    return (
        i_of_q,
        expected_k_vec,
        expected_k_error_vec,
        expected_corrected_intensities,
        expected_corrected_errors,
    )


def determine_common_mod_q_range(iqmod):
    """Determine the common Q1D range among all the wavelengths such that I(q, lambda) does exist

    Detailed requirement:
        Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process

    Parameters
    ----------
    iqmod :  ~drtsans.dataobjects.IQmod
        Input I(Q, wavelength) to find common Q range from

    Returns
    -------
    tuple
        minimum Q, maximum Q
    """
    # verify
    if iqmod.wavelength is None:
        raise RuntimeError("Input IQmod is required to have wavelength vector")

    # description of the algorithm for Q1D
    # get unique wave length
    wl_vec = np.unique(iqmod.wavelength)
    wl_vec.sort()

    # create a 2D array: wavelength, Q1D, intensity
    combo_vec = np.array([iqmod.wavelength, iqmod.mod_q, iqmod.intensity])
    # transpose for numpy selection and slicing
    combo_vec = combo_vec.transpose()

    # then for each unique wave length, determine the min and max
    min_q_vec = np.zeros_like(wl_vec, dtype=float)
    max_q_vec = np.zeros_like(wl_vec, dtype=float)
    for iw, wave_length in enumerate(wl_vec):
        # filter by wave length: column 0
        single_wl_q_vec = combo_vec[combo_vec[:, 0] == wave_length]
        # filter by intensity: column 2
        single_wl_q_vec = single_wl_q_vec[np.isfinite(single_wl_q_vec[:, 2])]
        # get minimum and max Q
        min_q_wl = np.min(single_wl_q_vec[:, 1])
        max_q_wl = np.max(single_wl_q_vec[:, 1])
        # set value
        min_q_vec[iw] = min_q_wl
        max_q_vec[iw] = max_q_wl

    # get the common qmin and qmax
    return np.max(min_q_vec), np.min(max_q_vec)


if __name__ == "__main__":
    pytest.main([__file__])
