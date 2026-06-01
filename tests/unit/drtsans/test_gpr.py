# standard imports
import os
import tempfile

# third party imports
import numpy as np
import pytest

# local imports
from drtsans.dataobjects import IQmod
from drtsans.extensions.gpr import (
    f_loglin,
    f_loglin_inv,
    f_loglin_deriv,
    run_gpr,
    create_png_plot,
    plotly_gpr_plot,
    generate_gpr_analysis,
    run_gpr_from_file,
)


class TestTransformFunctions:
    """Test log-linear transformation functions."""

    def test_f_loglin_basic(self):
        """Test basic log-linear transformation."""
        x = np.array([0.5, 1.0, 2.0, 4.0])
        result = f_loglin(x, x_c=1.0)
        # For x < x_c: (x - x_c) / x_c + log(x_c) = (x - 1) / 1 + 0 = x - 1
        # For x >= x_c: log(x)
        expected = np.array([-0.5, 0.0, np.log(2.0), np.log(4.0)])
        assert result == pytest.approx(expected, abs=1e-10)

    def test_f_loglin_inv_basic(self):
        """Test inverse log-linear transformation."""
        y = np.array([-0.5, 0.0, np.log(2.0), np.log(4.0)])
        result = f_loglin_inv(y, y_c=0.0)
        expected = np.array([0.5, 1.0, 2.0, 4.0])
        assert result == pytest.approx(expected, abs=1e-10)

    def test_f_loglin_roundtrip(self):
        """Test that f_loglin and f_loglin_inv are inverses."""
        x = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
        x_c = 1.0
        y_c = np.log(x_c)
        y = f_loglin(x, x_c=x_c)
        x_recovered = f_loglin_inv(y, y_c=y_c)
        assert x_recovered == pytest.approx(x, abs=1e-10)

    def test_f_loglin_deriv(self):
        """Test derivative of log-linear transformation."""
        x = np.array([0.5, 1.0, 2.0, 4.0])
        result = f_loglin_deriv(x, x_c=1.0)
        # For x < x_c: 1/x_c = 1
        # For x >= x_c: 1/x
        expected = np.array([1.0, 1.0, 0.5, 0.25])
        assert result == pytest.approx(expected, abs=1e-10)


class TestRunGPR:
    """Test core GPR analysis function."""

    @pytest.fixture
    def synthetic_data(self):
        """Create synthetic I(Q) data for testing."""
        # Generate Q values on log scale
        q = np.logspace(-2, 0, 50)  # 0.01 to 1.0

        # Create power-law I(Q) with noise
        I_true = 100 * q**-2.5
        noise_level = 0.1
        np.random.seed(42)
        noise = np.random.normal(0, noise_level * I_true, size=len(q))
        I = I_true + noise
        I_err = noise_level * I_true

        return q, I, I_err

    def test_run_gpr_basic(self, synthetic_data):
        """Test GPR with default parameters on synthetic data."""
        q, I, I_err = synthetic_data

        result = run_gpr(q, I, I_err)

        # Check that result has expected keys
        assert "q" in result
        assert "I_gpr" in result
        assert "I_gpr_err" in result
        assert "I_exp" in result
        assert "I_err_exp" in result
        assert "sigma_best" in result

        # Check array shapes match
        assert len(result["q"]) == len(q)
        assert len(result["I_gpr"]) == len(q)
        assert len(result["I_gpr_err"]) == len(q)

        # Check that GPR uncertainties are positive
        assert np.all(result["I_gpr_err"] > 0)

        # Check that sigma_best is reasonable
        assert 0.01 < result["sigma_best"] < 1000

    def test_run_gpr_with_log_transforms(self, synthetic_data):
        """Test GPR with log transformations enabled."""
        q, I, I_err = synthetic_data

        # Test with log Q (default)
        result_log_q = run_gpr(q, I, I_err, use_log_Q=True)
        assert "I_gpr" in result_log_q

        # Test without log Q
        result_no_log_q = run_gpr(q, I, I_err, use_log_Q=False)
        assert "I_gpr" in result_no_log_q

        # Test with log I
        result_log_i = run_gpr(q, I, I_err, use_log_I=True)
        assert "I_gpr" in result_log_i

    def test_run_gpr_with_dq(self, synthetic_data):
        """Test GPR with Q uncertainties provided."""
        q, I, I_err = synthetic_data
        dq = 0.01 * q  # 1% uncertainty in Q

        result = run_gpr(q, I, I_err, dq=dq)

        assert "dq" in result
        assert len(result["dq"]) == len(q)

    def test_run_gpr_invalid_input_lengths(self, synthetic_data):
        """Test that mismatched array lengths raise ValueError."""
        q, I, I_err = synthetic_data

        with pytest.raises(ValueError, match="must have the same length"):
            run_gpr(q[:-5], I, I_err)

        with pytest.raises(ValueError, match="must have the same length"):
            run_gpr(q, I[:-5], I_err)

    def test_run_gpr_negative_q(self, synthetic_data):
        """Test that negative Q values raise ValueError."""
        q, I, I_err = synthetic_data
        q_negative = q.copy()
        q_negative[0] = -0.01

        with pytest.raises(ValueError, match="Q values must be positive"):
            run_gpr(q_negative, I, I_err)

    def test_run_gpr_zero_errors_handled(self, synthetic_data):
        """Test that zero/negative errors are handled gracefully."""
        q, I, I_err = synthetic_data
        I_err_zero = I_err.copy()
        I_err_zero[5] = 0.0
        I_err_zero[10] = -0.1

        # Should not raise, but should log warning
        result = run_gpr(q, I, I_err_zero)
        assert "I_gpr" in result

    def test_run_gpr_custom_parameters(self, synthetic_data):
        """Test GPR with custom kernel and background parameters."""
        q, I, I_err = synthetic_data

        result = run_gpr(
            q,
            I,
            I_err,
            lmbda=0.5,  # larger kernel length scale
            background_filter_width=0.3,
            background_intensity_offset=2.0,
        )

        assert "I_gpr" in result
        assert len(result["I_gpr"]) == len(q)

    def test_run_gpr_small_dataset(self):
        """Test GPR with minimal dataset."""
        q = np.array([0.01, 0.1, 1.0])
        I = np.array([100.0, 10.0, 1.0])
        I_err = np.array([10.0, 1.0, 0.1])

        result = run_gpr(q, I, I_err)

        assert len(result["I_gpr"]) == 3
        assert np.all(np.isfinite(result["I_gpr"]))


class TestPlotting:
    """Test plotting functions."""

    @pytest.fixture
    def gpr_result(self):
        """Create mock GPR result for plotting tests."""
        q = np.logspace(-2, 0, 30)
        I_exp = 100 * q**-2.5
        I_err_exp = 0.1 * I_exp
        I_gpr = 100 * q**-2.5
        I_gpr_err = 0.05 * I_gpr

        return {
            "q": q,
            "I_gpr": I_gpr,
            "I_gpr_err": I_gpr_err,
            "I_exp": I_exp,
            "I_err_exp": I_err_exp,
            "sigma_best": 1.5,
        }

    def test_create_png_plot(self, gpr_result):
        """Test PNG plot creation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_gpr.png")

            result_path = create_png_plot(gpr_result, output_path, title="Test GPR")

            assert result_path == output_path
            assert os.path.isfile(output_path)
            assert os.path.getsize(output_path) > 1000  # Should be > 1KB

    def test_create_png_plot_no_title(self, gpr_result):
        """Test PNG plot creation without title."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_gpr_notitle.png")

            create_png_plot(gpr_result, output_path)

            assert os.path.isfile(output_path)

    def test_plotly_gpr_plot(self, gpr_result):
        """Test Plotly HTML plot generation."""
        html = plotly_gpr_plot(gpr_result, title="Test GPR", labels={"exp": "Experimental", "gpr": "GPR Fit"})

        assert isinstance(html, str)
        # Basic check that it's HTML-like
        if html:  # plot_publisher might not be available in all test environments
            assert "Q" in html or "div" in html.lower()

    def test_plotly_gpr_plot_default_labels(self, gpr_result):
        """Test Plotly plot with default labels."""
        html = plotly_gpr_plot(gpr_result)

        assert isinstance(html, str)


class TestGenerateGPRAnalysis:
    """Test entry point function for autoreduction integration."""

    @pytest.fixture
    def mock_iqmod(self):
        """Create mock IQmod object."""
        q = np.logspace(-2, 0, 30)
        I = 100 * q**-2.5
        I_err = 0.1 * I
        dq = 0.01 * q

        return IQmod(intensity=I, error=I_err, mod_q=q, delta_mod_q=dq)

    def test_generate_gpr_analysis_single_profile(self, mock_iqmod):
        """Test GPR analysis with single IQmod object."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html, png_files, dat_files = generate_gpr_analysis(mock_iqmod, tmpdir, "EQSANS_12345")

            assert isinstance(html, str)
            assert isinstance(png_files, list)
            assert isinstance(dat_files, list)

            # Check files were created
            assert len(png_files) == 1
            assert len(dat_files) == 1
            assert os.path.isfile(png_files[0])
            assert os.path.isfile(dat_files[0])

            # Check file naming
            assert "EQSANS_12345" in png_files[0]
            assert "EQSANS_12345" in dat_files[0]
            assert "_Iq_gpr.png" in png_files[0]
            assert "_Iq_gpr.dat" in dat_files[0]

    def test_generate_gpr_analysis_multiple_profiles(self, mock_iqmod):
        """Test GPR analysis with multiple IQmod objects (wedges)."""
        # Create 2 wedges
        iqmod_list = [mock_iqmod, mock_iqmod]

        with tempfile.TemporaryDirectory() as tmpdir:
            html, png_files, dat_files = generate_gpr_analysis(iqmod_list, tmpdir, "EQSANS_12345")

            # Should process both wedges
            assert len(png_files) == 2
            assert len(dat_files) == 2

            # Check wedge suffixes
            assert "_wedge_0" in png_files[0]
            assert "_wedge_1" in png_files[1]

    def test_generate_gpr_analysis_invalid_directory(self, mock_iqmod):
        """Test handling of invalid output directory."""
        # Use a path that can't be created (e.g., under /dev/null)
        invalid_dir = "/dev/null/nonexistent"

        html, png_files, dat_files = generate_gpr_analysis(mock_iqmod, invalid_dir, "EQSANS_12345")

        # Should return empty results on failure
        assert html == ""
        assert png_files == []
        assert dat_files == []

    def test_generate_gpr_analysis_creates_directory(self, mock_iqmod):
        """Test that output directory is created if it doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            new_dir = os.path.join(tmpdir, "subdir", "output")

            html, png_files, dat_files = generate_gpr_analysis(mock_iqmod, new_dir, "EQSANS_12345")

            assert os.path.isdir(new_dir)
            assert len(png_files) > 0


class TestRunGPRFromFile:
    """Test CLI file-based interface."""

    @pytest.fixture
    def test_dat_file(self):
        """Create test .dat file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("# Test I(Q) data\n")
            f.write("#Q (1/A)        I (1/cm)        dI (1/cm)       dQ (1/A)\n")
            q = np.logspace(-2, 0, 30)
            I = 100 * q**-2.5
            I_err = 0.1 * I
            dq = 0.01 * q
            for i in range(len(q)):
                f.write(f"{q[i]:.6E}\t{I[i]:.6E}\t{I_err[i]:.6E}\t{dq[i]:.6E}\n")
            filepath = f.name

        yield filepath

        # Cleanup
        if os.path.exists(filepath):
            os.remove(filepath)

    def test_run_gpr_from_file_basic(self, test_dat_file):
        """Test GPR from file with default output directory."""
        png_path, dat_path = run_gpr_from_file(test_dat_file)

        assert os.path.isfile(png_path)
        assert os.path.isfile(dat_path)
        assert png_path.endswith("_gpr.png")
        assert dat_path.endswith("_gpr.dat")

        # Cleanup
        os.remove(png_path)
        os.remove(dat_path)

    def test_run_gpr_from_file_custom_output_dir(self, test_dat_file):
        """Test GPR from file with custom output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            png_path, dat_path = run_gpr_from_file(test_dat_file, output_dir=tmpdir)

            assert tmpdir in png_path
            assert tmpdir in dat_path
            assert os.path.isfile(png_path)
            assert os.path.isfile(dat_path)

    def test_run_gpr_from_file_nonexistent(self):
        """Test that nonexistent file raises error."""
        with pytest.raises(FileNotFoundError):
            run_gpr_from_file("/nonexistent/file.dat")

    def test_run_gpr_from_file_invalid_format(self):
        """Test that invalid file format raises error."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            f.write("# Invalid data\n")
            f.write("# Only one column\n")
            f.write("1.0\n2.0\n3.0\n")
            filepath = f.name

        try:
            with pytest.raises(ValueError, match="must have at least 3 columns"):
                run_gpr_from_file(filepath)
        finally:
            os.remove(filepath)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
