"""Integration tests for GPR module with realistic workflows."""

# standard imports
import os
import tempfile

# third party imports
import numpy as np
import pytest

# local imports
from drtsans.dataobjects import IQmod
from drtsans.extensions.gpr import generate_gpr_analysis, run_gpr_from_file


@pytest.fixture
def realistic_iq_data():
    """Create realistic I(Q) data similar to EQSANS output.

    Uses a power-law decay with noise, typical of SANS data.
    """
    # Q range typical for EQSANS (0.005 to 0.5 Å^-1)
    q = np.logspace(np.log10(0.005), np.log10(0.5), 100)

    # Power-law: I(Q) = A * Q^(-alpha) with alpha ~ 4 (Porod scattering)
    # Plus a constant background
    A = 1e7
    alpha = 3.5
    background = 100
    I_true = A * q ** (-alpha) + background

    # Add realistic noise (5-10% relative error)
    relative_error = 0.08
    I_err = I_true * relative_error
    rng = np.random.default_rng(42)
    I = I_true + rng.normal(0, I_err)

    # Ensure no negative intensities
    I = np.maximum(I, I_err)

    # dQ typical for EQSANS (5% of Q)
    dq = 0.05 * q

    return IQmod(intensity=I, error=I_err, mod_q=q, delta_mod_q=dq)


@pytest.fixture
def temp_output_dir():
    """Create and cleanup temporary directory for test outputs."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    # Cleanup handled by system temp cleanup


class TestGPRWithRealisticData:
    """Test GPR analysis with realistic SANS-like data."""

    def test_single_profile_complete_workflow(self, realistic_iq_data, temp_output_dir):
        """Test complete GPR workflow with single I(Q) profile."""
        html_report, png_files, dat_files = generate_gpr_analysis(realistic_iq_data, temp_output_dir, "EQSANS_TEST")

        # Check HTML report generated
        assert isinstance(html_report, str)
        assert len(html_report) > 0
        assert "<div" in html_report  # Plotly div

        # Check files created
        assert len(png_files) == 1
        assert len(dat_files) == 1
        assert os.path.exists(png_files[0])
        assert os.path.exists(dat_files[0])

        # Verify PNG file is not empty
        assert os.path.getsize(png_files[0]) > 1000  # At least 1KB

        # Verify DAT file structure
        data = np.loadtxt(dat_files[0], skiprows=2)
        assert data.shape[1] == 4  # Q, I_GPR, I_GPR_err, dQ
        assert data.shape[0] > 0  # Has data points

        # Check Q values are monotonically increasing
        assert np.all(np.diff(data[:, 0]) > 0)

        # Check no NaN or Inf values
        assert not np.any(np.isnan(data))
        assert not np.any(np.isinf(data))

    def test_multiple_profiles_wedges(self, realistic_iq_data, temp_output_dir):
        """Test GPR with multiple profiles (simulating wedge analysis)."""
        # Create 4 wedge profiles with slightly different characteristics
        profiles = []
        for i in range(4):
            # Vary intensity slightly for each wedge
            scale_factor = 1.0 + 0.1 * (i - 1.5)
            wedge_profile = IQmod(
                intensity=realistic_iq_data.intensity * scale_factor,
                error=realistic_iq_data.error * scale_factor,
                mod_q=realistic_iq_data.mod_q,
                delta_mod_q=realistic_iq_data.delta_mod_q,
            )
            profiles.append(wedge_profile)

        html_report, png_files, dat_files = generate_gpr_analysis(profiles, temp_output_dir, "EQSANS_WEDGE_TEST")

        # Check multiple outputs created
        assert len(png_files) == 4
        assert len(dat_files) == 4

        # Verify all files exist
        for png_file, dat_file in zip(png_files, dat_files):
            assert os.path.exists(png_file)
            assert os.path.exists(dat_file)
            assert "_wedge_" in png_file
            assert "_wedge_" in dat_file

        # HTML should contain multiple plots
        assert html_report.count("<div") >= 4

    def test_gpr_parameters_effect(self, realistic_iq_data, temp_output_dir):
        """Test that GPR parameters affect the output."""
        # Run with default parameters
        _, _, dat_files_default = generate_gpr_analysis(
            realistic_iq_data, temp_output_dir, "EQSANS_DEFAULT", use_log_I=False, lmbda=0.25
        )

        # Run with different parameters
        _, _, dat_files_modified = generate_gpr_analysis(
            realistic_iq_data, temp_output_dir, "EQSANS_MODIFIED", use_log_I=True, lmbda=0.5
        )

        # Load both results
        data_default = np.loadtxt(dat_files_default[0], skiprows=2)
        data_modified = np.loadtxt(dat_files_modified[0], skiprows=2)

        # Results should differ (not identical)
        assert not np.allclose(data_default[:, 1], data_modified[:, 1], rtol=0.01)

    def test_low_statistics_data(self, temp_output_dir):
        """Test GPR with low statistics data (few points)."""
        # Only 10 data points
        q = np.logspace(np.log10(0.01), np.log10(0.3), 10)
        I = 1e6 * q ** (-3.5)
        I_err = I * 0.1

        iq_data = IQmod(intensity=I, error=I_err, mod_q=q, delta_mod_q=0.05 * q)

        html_report, png_files, dat_files = generate_gpr_analysis(iq_data, temp_output_dir, "EQSANS_LOWSTATS")

        # Should still produce output
        assert len(png_files) == 1
        assert len(dat_files) == 1
        assert os.path.exists(png_files[0])
        assert os.path.exists(dat_files[0])

    def test_high_noise_data(self, temp_output_dir):
        """Test GPR with high noise data."""
        q = np.logspace(np.log10(0.005), np.log10(0.5), 50)
        I_true = 1e6 * q ** (-3.5) + 100

        # 50% relative error (very noisy)
        I_err = I_true * 0.5
        I = I_true + np.random.normal(0, I_err)
        I = np.maximum(I, I_err * 0.1)  # Prevent negative intensities

        iq_data = IQmod(intensity=I, error=I_err, mod_q=q, delta_mod_q=0.05 * q)

        html_report, png_files, dat_files = generate_gpr_analysis(iq_data, temp_output_dir, "EQSANS_HIGHNOISE")

        # Should handle gracefully
        assert len(png_files) == 1
        assert len(dat_files) == 1

        # Verify output is reasonable
        data = np.loadtxt(dat_files[0], skiprows=2)
        assert not np.any(np.isnan(data))
        assert not np.any(np.isinf(data))


class TestGPRFileInterface:
    """Test GPR file-based interface."""

    def test_run_gpr_from_file_complete(self, realistic_iq_data, temp_output_dir):
        """Test run_gpr_from_file with realistic data file."""
        # Create input file
        input_file = os.path.join(temp_output_dir, "test_input_Iq.dat")
        header = "Q (1/A)\tI (1/cm)\tI_err (1/cm)\tdQ (1/A)"
        data = np.column_stack(
            (
                realistic_iq_data.mod_q,
                realistic_iq_data.intensity,
                realistic_iq_data.error,
                realistic_iq_data.delta_mod_q,
            )
        )
        np.savetxt(input_file, data, fmt="%.6e", delimiter="\t", header=header)

        # Run GPR from file
        png_path, dat_path = run_gpr_from_file(input_file, output_dir=temp_output_dir)

        # Check outputs
        assert os.path.exists(png_path)
        assert os.path.exists(dat_path)
        assert png_path.endswith("_Iq_gpr.png")
        assert dat_path.endswith("_Iq_gpr.dat")

        # Verify output data
        output_data = np.loadtxt(dat_path, skiprows=2)
        assert output_data.shape[1] == 4
        assert output_data.shape[0] > 0

    def test_run_gpr_from_file_custom_parameters(self, realistic_iq_data, temp_output_dir):
        """Test run_gpr_from_file with custom GPR parameters."""
        # Create input file
        input_file = os.path.join(temp_output_dir, "test_params_Iq.dat")
        header = "Q (1/A)\tI (1/cm)\tI_err (1/cm)\tdQ (1/A)"
        data = np.column_stack(
            (
                realistic_iq_data.mod_q,
                realistic_iq_data.intensity,
                realistic_iq_data.error,
                realistic_iq_data.delta_mod_q,
            )
        )
        np.savetxt(input_file, data, fmt="%.6e", delimiter="\t", header=header)

        # Run with custom parameters
        png_path, dat_path = run_gpr_from_file(
            input_file, output_dir=temp_output_dir, use_log_I=True, use_log_Q=True, lmbda=0.5
        )

        assert os.path.exists(png_path)
        assert os.path.exists(dat_path)


class TestGPRErrorHandling:
    """Test GPR error handling in integration scenarios."""

    def test_empty_profile_list(self, temp_output_dir):
        """Test handling of empty profile list."""
        # Empty list should be handled gracefully, not crash
        html_report, png_files, dat_files = generate_gpr_analysis([], temp_output_dir, "EQSANS_EMPTY")

        # Should return empty results without crashing
        assert html_report == ""
        assert len(png_files) == 0
        assert len(dat_files) == 0

    def test_invalid_output_directory(self, realistic_iq_data):
        """Test handling of invalid output directory."""
        invalid_dir = "/invalid/nonexistent/directory/path"

        # Should create directory or handle gracefully
        html_report, png_files, dat_files = generate_gpr_analysis(realistic_iq_data, invalid_dir, "EQSANS_INVALID")

        # Should return empty results on failure
        assert html_report == ""
        assert len(png_files) == 0
        assert len(dat_files) == 0

    def test_mismatched_data_lengths(self, temp_output_dir):
        """Test handling of mismatched Q and I arrays."""
        # IQmod validates dimensions during construction
        with pytest.raises(TypeError, match="Shape mismatch"):
            IQmod(
                intensity=np.array([1, 2, 3]),
                error=np.array([0.1, 0.2, 0.3]),
                mod_q=np.array([0.01, 0.02]),  # Wrong length
                delta_mod_q=np.array([0.001, 0.002]),
            )
