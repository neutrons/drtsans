# Test drtsans.tof.eqsans.incoherence_correction_1d
import pytest


def test_incoherence_inelastic_correction():

    # generate testing data
    test_iq1d = generate_test_data()



def generate_test_data():
    # Generate test data given in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/
    # /b3b4038f44443385afe4252bb2316be3/inelastic_incoherent_avg_example.xlsx
    # Denoted as TEST1

    # Intensity vector
    intensity_vec = np.array([
        0.1, 0.13, np.nan, np.nan, np.nan,  # q = 0.01
        0.1, 0.13, 0.15, np.nan, np.nan,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        0.1, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, 0.13, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, 0.15, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])

    # Q vector
    vec_q = np.arange(1, 21) * 0.01
    vec_q = np.repeat(vec_q, 5)
    print(f'q: {vec_q.shape}')

    # Wavelength vector
    wavelength_vec = np.arange(1, 5) * 1.
    wavelength_vec = np.tile(wavelength_vec, 20)
    print(f'lambda: {wavelength_vec.shape}')

    # Error
    error_vec = np.sqrt(intensity_vec)

    # Construct IQmod
    i_of_q = IQmod(intensity=intensity_vec,
                   error=error_vec,
                   mod_q=vec_q,
                   wavelength=wavelength_vec)

    return i_of_q


def generate_expected_b_factors:
    # Generate test data given in
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/
    # /b3b4038f44443385afe4252bb2316be3/inelastic_incoherent_avg_example.xlsx
    # Denoted as TEST1

    # Expected B vectors
    b_factor_vec = np.array([np.nan, -0.03, -0.05, -0.04, -0.01])

    return b_factor_vec


def generate_expected_corrected_intensities:

    # Expected corrected intensities
    corrected_intensity_vec = np.array([
        0.1, 0.1, np.nan, np.nan, np.nan,  # q = 0.01
        0.1, 0.1, 0.1, np.nan, np.nan,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1, 0.1,
        np.nan, 0.1, 0.1, 0.1, 0.1,
        np.nan, 0.1, 0.1, 0.1, 0.1,
        np.nan, np.nan, 0.1, 0.1, 0.1,
        np.nan, np.nan, 0.1, 0.1, 0.1,
        np.nan, np.nan, np.nan, 0.1, 0.1,
        np.nan, np.nan, np.nan, 0.14, 0.11,
        np.nan, np.nan, np.nan, np.nan, 0.11,
    ])

    return corrected_intensity_vec



if __name__ == '__main__':
    pytest.main([__file__])

