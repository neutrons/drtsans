"""Integration test for EWM-13940: Scalar and 360° wedge binning consistency with inelastic correction

This test verifies that scalar I(Q) and 360° wedge I(Q) produce identical results when
inelastic/incoherent correction is enabled.

EWM-13940: Inelastic/incoherent correction produces different results for scalar and wedge binning
Fix: Apply corrections to unbinned data before mode-specific binning
"""

import os

import numpy as np
import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import DeleteWorkspace, mtd

from drtsans.tof.eqsans import reduction_parameters, load_all_files, reduce_single_configuration


@pytest.mark.datarepo
def test_scalar_vs_wedge360_with_inelastic_correction(datarepo_dir, temp_directory):
    """Test that scalar and symmetric wedge I(Q) produce consistent results with inelastic correction.

    This is the core regression test for EWM-13940. Before the fix, scalar and wedge binning
    produced different results when inelastic correction was enabled because corrections were
    applied to intermediate binned data. After the fix, corrections are applied to unbinned
    data before binning, ensuring consistent results.

    """

    # Base configuration using existing test data
    base_configuration = {
        "instrumentName": "EQSANS",
        "iptsNumber": "27799",
        "sample": {"runNumber": "125707", "thickness": 1, "transmission": {"runNumber": "", "value": "1"}},
        "background": {"runNumber": "", "transmission": {"runNumber": "", "value": ""}},
        "emptyTransmission": {"runNumber": "125701", "value": ""},
        "beamCenter": {"runNumber": "125701"},
        "configuration": {
            "cutTOFmax": 2000.0,
            "cutTOFmin": 500.0,
            "wavelengthStep": 0.1,
            "wavelengthStepType": "constant Delta lambda",
            "sampleApertureSize": 10,
            "numQBins": 100,
            "numQxQyBins": 80,
            "QbinType": "linear",
            "Qmin": None,
            "Qmax": None,
            # Enable inelastic correction
            "fitInelasticIncoh": [True],
            "selectMinIncoh": True,
            # Enable elastic correction for completeness
            "elasticReference": {
                "runNumber": "124680",
                "thickness": "1.0",
                "transmission": {"runNumber": None, "value": "1.0"},
            },
            "elasticReferenceBkgd": {
                "runNumber": None,
                "transmission": {"runNumber": None, "value": "0.9"},
            },
            # Files
            "darkFileName": "/bin/true",
            "useDefaultMask": True,
            "normalization": "Total charge",
            "detectorOffset": "80",
            "sampleOffset": "314.5",
        },
    }

    # Set up data directories
    datadir = os.path.join(datarepo_dir.eqsans, "test_corrections")
    base_configuration["dataDirectories"] = datadir
    base_configuration["configuration"]["maskFileName"] = os.path.join(datadir, "beamstop_mask_4m_ext.nxs")
    base_configuration["configuration"]["sensitivityFileName"] = os.path.join(
        datadir, "Sensitivity_patched_thinPMMA_4m_124972.nxs"
    )
    base_configuration["configuration"]["instrumentConfigurationDir"] = os.path.join(
        datarepo_dir.eqsans, "instrument_configuration"
    )
    base_configuration["configuration"]["beamFluxFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt"
    )

    # Create temp output directory
    test_dir = temp_directory(prefix="test_ewm13940")

    # Test 1: Scalar binning
    print("\n=== Running reduction with SCALAR binning ===")
    config_scalar = base_configuration.copy()
    config_scalar["configuration"] = base_configuration["configuration"].copy()
    config_scalar["configuration"]["outputDir"] = test_dir
    config_scalar["outputFileName"] = "EWM13940_scalar"
    config_scalar["configuration"]["1DQbinType"] = "scalar"

    with amend_config(data_dir=datarepo_dir.eqsans):
        input_config_scalar = reduction_parameters(config_scalar)
        input_config_scalar["configuration"]["darkFileName"] = None
        loaded_scalar = load_all_files(input_config_scalar)
        output_scalar = reduce_single_configuration(loaded_scalar, input_config_scalar)

    assert output_scalar is not None
    assert len(output_scalar) == 1
    iq_scalar = output_scalar[0].I1D_main[0]

    # Test 2: Symmetric 179° wedge binning
    # Single wedge 0-179° with autoSymmetricWedges=True covers full azimuthal range
    # Note: Wedge angle must be <180° (not <=180°)
    print("\n=== Running reduction with symmetric 179° WEDGE binning ===")
    config_wedge = base_configuration.copy()
    config_wedge["configuration"] = base_configuration["configuration"].copy()
    config_wedge["configuration"]["outputDir"] = test_dir
    config_wedge["outputFileName"] = "EWM13940_wedge_symmetric"
    config_wedge["configuration"]["1DQbinType"] = "wedge"
    config_wedge["configuration"]["WedgeMinAngles"] = [0]
    config_wedge["configuration"]["WedgeMaxAngles"] = [179]
    config_wedge["configuration"]["autoSymmetricWedges"] = True

    with amend_config(data_dir=datarepo_dir.eqsans):
        input_config_wedge = reduction_parameters(config_wedge)
        input_config_wedge["configuration"]["darkFileName"] = None
        loaded_wedge = load_all_files(input_config_wedge)
        output_wedge = reduce_single_configuration(loaded_wedge, input_config_wedge)

    assert output_wedge is not None
    assert len(output_wedge) == 1
    # For symmetric wedge, we get one I(Q) that combines data from symmetric 179° wedges
    iq_wedge_symmetric = output_wedge[0].I1D_main[0]

    # Compare scalar vs symmetric 179° wedge
    print("\n=== Comparing scalar vs symmetric 179° wedge results ===")

    # Check that Q bins are identical
    np.testing.assert_allclose(
        iq_scalar.mod_q,
        iq_wedge_symmetric.mod_q,
        rtol=1e-10,
        err_msg="Q bins should be identical for scalar and symmetric 179° wedge",
    )

    # Check that intensities match
    # Allow small numerical differences from floating point operations
    intensity_diff = np.abs(iq_scalar.intensity - iq_wedge_symmetric.intensity)
    intensity_mean = np.abs(iq_scalar.intensity)

    # Calculate relative differences (avoiding division by near-zero)
    mask_nonzero = intensity_mean > 1e-10
    rel_diff = np.zeros_like(intensity_diff)
    rel_diff[mask_nonzero] = intensity_diff[mask_nonzero] / intensity_mean[mask_nonzero]

    mean_rel_diff = np.mean(rel_diff[mask_nonzero])
    max_rel_diff = np.max(rel_diff[mask_nonzero])

    print(f"Mean relative difference: {mean_rel_diff:.6e}")
    print(f"Max relative difference: {max_rel_diff:.6e}")

    # Assert that differences are within acceptable tolerance
    # Scalar and symmetric wedge aren't perfectly identical due to:
    # - Wedge covers 179° x 2 = 358° (missing 2° vs scalar's 360°)
    # - Different binning implementations
    # But they should be very close after the EWM-13940 fix (<0.5% mean, <1% max)
    assert mean_rel_diff < 5e-3, (
        f"Mean relative difference {mean_rel_diff:.2e} exceeds 0.5% threshold. "
        f"Scalar and symmetric 179° wedge I(Q) should be nearly identical after EWM-13940 fix."
    )
    assert max_rel_diff < 1e-2, (
        f"Max relative difference {max_rel_diff:.2e} exceeds 1% threshold. "
        f"Scalar and symmetric 179° wedge I(Q) should be nearly identical after EWM-13940 fix."
    )

    # Also check errors match
    error_diff = np.abs(iq_scalar.error - iq_wedge_symmetric.error)
    error_mean = np.abs(iq_scalar.error)
    mask_error_nonzero = error_mean > 1e-10
    rel_error_diff = np.zeros_like(error_diff)
    rel_error_diff[mask_error_nonzero] = error_diff[mask_error_nonzero] / error_mean[mask_error_nonzero]

    mean_rel_error_diff = np.mean(rel_error_diff[mask_error_nonzero])
    max_rel_error_diff = np.max(rel_error_diff[mask_error_nonzero])

    print(f"Mean relative error difference: {mean_rel_error_diff:.6e}")
    print(f"Max relative error difference: {max_rel_error_diff:.6e}")

    assert mean_rel_error_diff < 5e-3, (
        f"Mean relative error difference {mean_rel_error_diff:.2e} exceeds 0.5% threshold"
    )
    assert max_rel_error_diff < 1e-2, f"Max relative error difference {max_rel_error_diff:.2e} exceeds 1% threshold"

    print("\n✓ SUCCESS: Scalar and symmetric 179° wedge I(Q) are nearly identical (EWM-13940 fix validated)")

    # Cleanup
    cleanup_workspaces()


@pytest.mark.datarepo
def test_multiple_wedges_with_inelastic_correction(datarepo_dir, temp_directory):
    """Test that multiple wedges can be binned consistently with inelastic correction.

    This test verifies that the fix doesn't break multi-wedge binning scenarios.
    """

    # Base configuration
    configuration = {
        "instrumentName": "EQSANS",
        "iptsNumber": "27799",
        "sample": {"runNumber": "125707", "thickness": 1, "transmission": {"runNumber": "", "value": "1"}},
        "background": {"runNumber": "", "transmission": {"runNumber": "", "value": ""}},
        "emptyTransmission": {"runNumber": "125701", "value": ""},
        "beamCenter": {"runNumber": "125701"},
        "configuration": {
            "cutTOFmax": 2000.0,
            "cutTOFmin": 500.0,
            "wavelengthStep": 0.1,
            "wavelengthStepType": "constant Delta lambda",
            "sampleApertureSize": 10,
            "numQBins": 100,
            "numQxQyBins": 80,
            "1DQbinType": "wedge",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "QbinType": "linear",
            # Enable inelastic correction
            "fitInelasticIncoh": [True],
            "selectMinIncoh": True,
            "elasticReference": {
                "runNumber": "124680",
                "thickness": "1.0",
                "transmission": {"runNumber": None, "value": "1.0"},
            },
            "elasticReferenceBkgd": {
                "runNumber": None,
                "transmission": {"runNumber": None, "value": "0.9"},
            },
            "darkFileName": "/bin/true",
            "useDefaultMask": True,
            "normalization": "Total charge",
            "detectorOffset": "80",
            "sampleOffset": "314.5",
        },
    }

    # Set up data directories
    datadir = os.path.join(datarepo_dir.eqsans, "test_corrections")
    configuration["dataDirectories"] = datadir
    configuration["configuration"]["maskFileName"] = os.path.join(datadir, "beamstop_mask_4m_ext.nxs")
    configuration["configuration"]["sensitivityFileName"] = os.path.join(
        datadir, "Sensitivity_patched_thinPMMA_4m_124972.nxs"
    )
    configuration["configuration"]["instrumentConfigurationDir"] = os.path.join(
        datarepo_dir.eqsans, "instrument_configuration"
    )
    configuration["configuration"]["beamFluxFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt"
    )

    test_dir = temp_directory(prefix="test_ewm13940_wedges")
    configuration["configuration"]["outputDir"] = test_dir
    configuration["outputFileName"] = "EWM13940_multi_wedge"

    with amend_config(data_dir=datarepo_dir.eqsans):
        input_config = reduction_parameters(configuration)
        input_config["configuration"]["darkFileName"] = None
        loaded = load_all_files(input_config)
        output = reduce_single_configuration(loaded, input_config)

    assert output is not None
    assert len(output) == 1
    # Should have 2 wedges
    assert len(output[0].I1D_main) == 2

    print(f"✓ Multi-wedge binning produced {len(output[0].I1D_main)} wedges as expected")

    # Cleanup
    cleanup_workspaces()


def cleanup_workspaces():
    """Clean up Mantid workspaces created during tests"""
    workspaces_to_delete = []
    for ws_name in mtd.getObjectNames():
        if any(
            ws_name.startswith(prefix)
            for prefix in [
                "_empty",
                "_mask",
                "_sensitivity",
                "processed_",
                "_EQSANS_",
                "EQ_TEST_",
                "_bkgd_",
                "_sample_",
            ]
        ):
            workspaces_to_delete.append(ws_name)

    for ws_name in workspaces_to_delete:
        if mtd.doesExist(ws_name):
            try:
                DeleteWorkspace(ws_name)
            except Exception as e:
                print(f"Warning: Could not delete workspace {ws_name}: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
