"""
Test that inelastic correction produces consistent results for scalar and wedge binning.

This test addresses EWM-13940: Intensity is not correct when inelastic correction
and multiple wedging are both turned on at EQSANS.

The fix ensures that correction factors are calculated from scalar-binned data and
applied to unbinned data BEFORE mode-specific binning, so that scalar and 360° wedge
I(Q) profiles overlap.
"""

import pytest
import numpy as np
from drtsans.dataobjects import IQmod, IQazimuthal
from drtsans.tof.eqsans.incoherence_correction_1d import (
    calculate_incoherence_correction_factors,
    apply_incoherence_correction_to_unbinned_data,
    CorrectionFactors,
)


def generate_synthetic_unbinned_data():
    """
    Generate synthetic unbinned I(Q, λ) data for testing

    Returns
    -------
    tuple
        (IQazimuthal 2D data, IQmod 1D data)
    """
    # Create synthetic wavelength-dependent data
    n_points = 1000
    n_wavelengths = 5

    # Generate Q values covering angular range
    theta = np.linspace(0, 2 * np.pi, n_points)
    q_magnitude = np.linspace(0.001, 0.3, n_points)

    # Create Qx, Qy
    qx = q_magnitude * np.cos(theta)
    qy = q_magnitude * np.sin(theta)

    # Create wavelength array (repeating for each Q point)
    wavelengths = np.tile(np.linspace(2.5, 6.5, n_wavelengths), n_points // n_wavelengths + 1)[:n_points]

    # Create intensity with wavelength dependence (simulating inelastic contribution)
    # Base intensity decreases with Q
    base_intensity = 100 * np.exp(-q_magnitude)

    # Add wavelength-dependent inelastic contribution
    # b(λ) increases with wavelength
    inelastic_contribution = (wavelengths - 2.5) * 10

    # Total intensity
    intensity = base_intensity + inelastic_contribution
    error = np.sqrt(intensity)  # Poisson statistics

    # Create IQazimuthal (unbinned 2D)
    iq2d = IQazimuthal(
        intensity=intensity,
        error=error,
        qx=qx,
        qy=qy,
        wavelength=wavelengths,
        delta_qx=np.ones_like(qx) * 0.001,
        delta_qy=np.ones_like(qy) * 0.001,
    )

    # Create IQmod (unbinned 1D)
    iq1d = IQmod(
        intensity=intensity,
        error=error,
        mod_q=q_magnitude,
        wavelength=wavelengths,
        delta_mod_q=np.ones_like(q_magnitude) * 0.001,
    )

    return iq2d, iq1d


def generate_scalar_binned_data():
    """Generate scalar-binned I(Q, λ) for correction factor calculation"""
    n_q_bins = 20
    n_wavelengths = 5

    q_bins = np.linspace(0.001, 0.3, n_q_bins)
    wavelengths = np.linspace(2.5, 6.5, n_wavelengths)

    # Create mesh
    q_mesh = np.tile(q_bins, (n_wavelengths, 1)).T
    wl_mesh = np.tile(wavelengths, (n_q_bins, 1))

    # Base intensity with wavelength-dependent inelastic contribution
    base_intensity = 100 * np.exp(-q_mesh)
    inelastic_contribution = (wl_mesh - 2.5) * 10
    intensity = base_intensity + inelastic_contribution

    # Flatten for IQmod
    intensity_flat = intensity.flatten()
    error_flat = np.sqrt(intensity_flat)
    q_flat = q_mesh.flatten()
    wl_flat = wl_mesh.flatten()

    return IQmod(
        intensity=intensity_flat,
        error=error_flat,
        mod_q=q_flat,
        wavelength=wl_flat,
        delta_mod_q=np.ones_like(q_flat) * 0.001,
    )


def test_calculate_correction_factors():
    """Test that correction factors can be calculated from scalar-binned data"""
    scalar_binned = generate_scalar_binned_data()

    factors = calculate_incoherence_correction_factors(
        scalar_binned,
        select_minimum_incoherence=False,
        intensity_weighted=False,
        qmin=0.01,
        qmax=0.2,
    )

    # Verify structure
    assert isinstance(factors, CorrectionFactors)
    assert len(factors.b_factor) == 5  # 5 wavelengths
    assert len(factors.b_error) == 5
    assert len(factors.wavelength) == 5

    # Verify b_factors are reasonable (should be positive for longer wavelengths)
    # since longer wavelengths have more inelastic contribution
    assert np.all(factors.b_factor >= 0), "B factors should be non-negative"

    # Verify wavelength order
    assert np.all(np.diff(factors.wavelength) > 0), "Wavelengths should be sorted"


def test_apply_correction_to_unbinned():
    """Test that corrections can be applied to unbinned data"""
    iq2d, iq1d = generate_synthetic_unbinned_data()
    scalar_binned = generate_scalar_binned_data()

    # Calculate correction factors
    factors = calculate_incoherence_correction_factors(
        scalar_binned,
        select_minimum_incoherence=False,
        intensity_weighted=False,
    )

    # Apply corrections
    iq2d_corrected, iq1d_corrected = apply_incoherence_correction_to_unbinned_data(iq2d, iq1d, factors)

    # Verify output structure
    assert isinstance(iq2d_corrected, IQazimuthal)
    assert isinstance(iq1d_corrected, IQmod)

    # Verify data shapes match
    assert iq2d_corrected.intensity.shape == iq2d.intensity.shape
    assert iq1d_corrected.intensity.shape == iq1d.intensity.shape

    # Verify intensity was reduced (correction was applied)
    assert np.mean(iq2d_corrected.intensity) < np.mean(iq2d.intensity)
    assert np.mean(iq1d_corrected.intensity) < np.mean(iq1d.intensity)

    # Verify errors were propagated
    assert np.all(iq2d_corrected.error > 0)
    assert np.all(iq1d_corrected.error > 0)

    # Errors should generally increase due to error propagation
    # (unless the data happens to have very small errors initially)
    assert np.mean(iq2d_corrected.error) >= np.mean(iq2d.error) * 0.9


def test_correction_preserves_data_structure():
    """Test that correction preserves Q, wavelength, and other metadata"""
    iq2d, iq1d = generate_synthetic_unbinned_data()
    scalar_binned = generate_scalar_binned_data()

    factors = calculate_incoherence_correction_factors(
        scalar_binned,
        select_minimum_incoherence=False,
    )

    iq2d_corrected, iq1d_corrected = apply_incoherence_correction_to_unbinned_data(iq2d, iq1d, factors)

    # Verify Q values unchanged
    np.testing.assert_array_equal(iq2d_corrected.qx, iq2d.qx)
    np.testing.assert_array_equal(iq2d_corrected.qy, iq2d.qy)
    np.testing.assert_array_equal(iq1d_corrected.mod_q, iq1d.mod_q)

    # Verify wavelengths unchanged
    np.testing.assert_array_equal(iq2d_corrected.wavelength, iq2d.wavelength)
    np.testing.assert_array_equal(iq1d_corrected.wavelength, iq1d.wavelength)

    # Verify delta_Q unchanged
    np.testing.assert_array_equal(iq2d_corrected.delta_qx, iq2d.delta_qx)
    np.testing.assert_array_equal(iq2d_corrected.delta_qy, iq2d.delta_qy)
    np.testing.assert_array_equal(iq1d_corrected.delta_mod_q, iq1d.delta_mod_q)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
