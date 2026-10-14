"""Manual EQSANS sensitivity validation using the IPTS-37618 water run.

Run this test explicitly with::

    pixi run pytest -vv -m mount_eqsans \
        tests/manual/test_prepare_sensitivities_eqsans_water.py

The test requires access to ``/SNS/EQSANS/IPTS-37618/nexus`` and validates the
scientific result of the sensitivity map rather than comparing NeXus files
byte-for-byte.  It generates a sensitivity map, uses it in the H2O reduction,
and checks the high-Q flatness of every wavelength-band profile before the
inelastic incoherent correction.
"""

from pathlib import Path

import numpy as np
import pytest
from mantid.api import mtd
from mantid.kernel import amend_config

from drtsans.path import load_module
from drtsans.redparams import reduction_parameters
from drtsans.tof.eqsans import load_all_files, reduce_single_configuration


DATA_DIRECTORY = Path("/SNS/EQSANS/IPTS-37618/nexus")
INSTRUMENT_CONFIGURATION_DIRECTORY = Path("/SNS/EQSANS/shared/instrument_configuration")
SENSITIVITY_TOML = """
FLOOD_RUNS = [188982]
DIRECT_BEAM_RUNS = [188969]

BLOCKED_BEAM_RUNS = [188984]
EMPTY_CELL_RUNS = [188978]

TRANSMISSION_FLOOD_RUNS = [188974]
TRANSMISSION_EMPTY_CELL_RUNS = [188970]
TRANSMISSION_REFERENCE_RUNS = [188969]

SAMPLE_OFFSET = 285.0
DETECTOR_OFFSET = 66.714
LOW_TOF_CLIP = 1650.0
HIGH_TOF_CLIP = 3150.0

MASKED_PIXELS = "1-15,242-256"
MASK_BEAM_CENTER_RADIUS = 40.0
PIXEL_CALIBRATION = false
SCALE_COMPONENTS = { detector1 = [1.004124, 1.057996, 1.0] }
SOLID_ANGLE_CORRECTION = true

MOVING_DETECTORS = false
MIN_THRESHOLD = 0.1
MAX_THRESHOLD = 3.0
"""

MAX_HIGH_Q_RMS = 0.02
MAX_HIGH_Q_TREND = 0.02
HIGH_Q_POINTS = 10
MAX_RELATIVE_ERROR = 0.02


def _sensitivity_script():
    """Load the source sensitivity script used by the manual test."""
    repository_root = Path(__file__).parents[2]
    return load_module(repository_root / "scripts/prepare_sensitivities_eqsans.py")


def _reduction_input(sensitivity_file: Path, output_directory: Path) -> dict:
    """Create the H2O reduction configuration for the generated map."""
    return {
        "instrumentName": "EQSANS",
        "iptsNumber": "37618",
        "outputFileName": "H2O_1p3m_manual_sensitivity",
        "sample": {
            "runNumber": "188982",
            "thickness": 0.1,
            "transmission": {"runNumber": "188974", "value": ""},
        },
        "background": {
            "runNumber": "188978",
            "transmission": {"runNumber": "", "value": "1.0"},
        },
        "emptyTransmission": {"runNumber": "188970", "value": None},
        "beamCenter": {"runNumber": "188970"},
        "dataDirectories": str(DATA_DIRECTORY),
        "configuration": {
            "outputDir": str(output_directory),
            "instrumentConfigurationDir": str(INSTRUMENT_CONFIGURATION_DIRECTORY),
            "cutTOFmin": 1650.0,
            "cutTOFmax": 3150.0,
            "wavelengthStep": 0.1,
            "wavelengthStepType": "constant Delta lambda",
            "sampleOffset": 285.0,
            "useDetectorOffset": True,
            "detectorOffset": 66.714,
            "sampleApertureSize": 10,
            "maskFileName": "/SNS/EQSANS/shared/script/eqsanstools/mask_4m.nxs",
            "useDefaultMask": True,
            "darkFileName": "/bin/true",
            "normalization": "Total charge",
            "blockedBeamRunNumber": "188984",
            "beamFluxFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2026B_mp/bl6_flux_2026B_aug_rebinned.txt",
            "sensitivityFileName": str(sensitivity_file),
            "useSolidAngleCorrection": True,
            "useThetaDepTransCorrection": True,
            "mmRadiusForTransmission": 25,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 1.0,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numQBins": 100,
            "AnnularAngleBin": 5.0,
            "useErrorWeighting": True,
            "fitInelasticIncoh": True,
            "selectMinIncoh": True,
            "outputWavelengthDependentProfile": True,
            "incohfit_intensityweighted": False,
            "elasticReference": {
                "runNumber": "",
                "thickness": "1.0",
                "transmission": {"runNumber": None, "value": "0.9"},
            },
            "elasticReferenceBkgd": {
                "runNumber": "",
                "transmission": {"runNumber": None, "value": "0.9"},
            },
            "scaleComponents": {"detector1": [1.004124, 1.057996, 1.0]},
        },
    }


def _profile_flatness(profile_file: Path) -> tuple[float, float]:
    """Return weighted high-Q RMS and end-to-end trend for one profile."""
    profile = np.loadtxt(profile_file, comments="#", ndmin=2)
    if profile.shape[1] < 3:
        raise AssertionError(f"Profile has fewer than three columns: {profile_file}")

    high_q = profile[np.argsort(profile[:, 0])[-HIGH_Q_POINTS:], :3]
    q, intensity, uncertainty = high_q.T
    valid = (
        np.isfinite(q)
        & np.isfinite(intensity)
        & np.isfinite(uncertainty)
        & (intensity > 0.0)
        & (uncertainty > 0.0)
        & (uncertainty / intensity <= MAX_RELATIVE_ERROR)
    )
    if np.count_nonzero(valid) < 3:
        raise AssertionError(f"Insufficient reliable high-Q points in {profile_file}")

    q = q[valid]
    intensity = intensity[valid]
    uncertainty = uncertainty[valid]
    weights = 1.0 / uncertainty**2
    weighted_mean = np.average(intensity, weights=weights)
    rms = np.sqrt(np.average(((intensity - weighted_mean) / weighted_mean) ** 2, weights=weights))

    slope, intercept = np.polyfit(q, intensity, 1, w=1.0 / uncertainty)
    fitted_endpoints = slope * np.array([q.min(), q.max()]) + intercept
    trend = abs(fitted_endpoints[1] - fitted_endpoints[0]) / weighted_mean
    return float(rms), float(trend)


@pytest.mark.mount_eqsans
def test_eqsans_sensitivity_produces_flat_h2o_profiles(has_sns_mount, tmp_path):
    """Check that the generated EQSANS map produces flat high-Q water profiles."""
    # This is deliberately a manual test: the experiment data live on the SNS mount.
    if not has_sns_mount or not DATA_DIRECTORY.is_dir():
        pytest.skip("SNS EQSANS data mount is not available")

    # Check all runs before starting Mantid so a partial reduction is not mistaken for a test result.
    required_runs = (188969, 188970, 188974, 188978, 188982, 188984)
    missing_runs = [run for run in required_runs if not (DATA_DIRECTORY / f"EQSANS_{run}.nxs.h5").is_file()]
    if missing_runs:
        pytest.skip(f"Required IPTS-37618 runs are not available: {missing_runs}")

    # Keep the reviewer-supplied TOML input inside the test while sending all generated files to pytest's temp area.
    sensitivity_script = _sensitivity_script()
    toml_file = tmp_path / "eqsans_sensitivity.toml"
    sensitivity_file = tmp_path / "Sensitivity_patched_H2O_1o3m_188982.nxs"
    reduction_directory = tmp_path / "reduction"
    toml_file.write_text(SENSITIVITY_TOML)

    try:
        with amend_config(data_dir=str(DATA_DIRECTORY)):
            # First create the sensitivity map from the same flood, beam-center, background,
            # and transmission runs used by the manual validation recipe.
            sensitivity_script.main([str(toml_file), str(sensitivity_file)])
            assert sensitivity_file.is_file()

            # Reduce water with the generated map.  The reduction writes one I(Q) profile per
            # wavelength band; dark-current subtraction stays disabled to isolate sensitivity effects.
            reduction_input = reduction_parameters(
                _reduction_input(sensitivity_file, reduction_directory),
                validate=True,
            )
            # The path above satisfies schema validation while retaining the
            # requested manual-test condition: no dark-current subtraction.
            reduction_input["configuration"]["darkFileName"] = None
            loaded = load_all_files(reduction_input, prefix="MANUAL_EQSANS_SENSITIVITY")
            reduce_single_configuration(loaded, reduction_input)

        # The profiles before the inelastic incoherent correction measure the effect of
        # sensitivity normalization without conflating it with the later water correction.
        profiles = sorted(reduction_directory.rglob("IQ_*_before_b_correction.dat"))
        assert len(profiles) > 1, "The reduction did not produce wavelength-band I(Q) profiles"

        # Evaluate every wavelength band in the high-Q tail against the flat-water criterion.
        flatness = [_profile_flatness(profile) for profile in profiles]
        rms_values = [rms for rms, _ in flatness]
        trend_values = [trend for _, trend in flatness]
        assert max(rms_values) <= MAX_HIGH_Q_RMS, (
            f"Maximum weighted high-Q RMS is {max(rms_values):.3%}; limit is {MAX_HIGH_Q_RMS:.3%}"
        )
        assert max(trend_values) <= MAX_HIGH_Q_TREND, (
            f"Maximum weighted high-Q trend is {max(trend_values):.3%}; limit is {MAX_HIGH_Q_TREND:.3%}"
        )
    finally:
        # The reduction creates many Mantid workspaces; clear them before pytest continues.
        mtd.clear()
