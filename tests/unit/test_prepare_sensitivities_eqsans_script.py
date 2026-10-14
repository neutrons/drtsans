from pathlib import Path

import pytest

from drtsans.path import load_module


_ROOT_DIR = Path(__file__).parents[2]
_SCRIPT = load_module(_ROOT_DIR / "scripts/prepare_sensitivities_eqsans.py")


def test_configuration_accepts_scalar_and_sequence_run_inputs():
    configuration = _SCRIPT.merge_configuration(
        {
            "FLOOD_RUNS": 1234,
            "DIRECT_BEAM_RUNS": ["1235"],
        }
    )

    assert configuration["FLOOD_RUNS"] == (1234,)
    assert configuration["DIRECT_BEAM_RUNS"] == (1235,)


def test_configuration_preserves_named_runs_and_paths():
    configuration = _SCRIPT.merge_configuration(
        {
            "FLOOD_RUNS": ["EQSANS_1234", "/path/to/EQSANS_1235.nxs.h5"],
            "DIRECT_BEAM_RUNS": "EQSANS_1236",
        }
    )

    assert configuration["FLOOD_RUNS"] == ("EQSANS_1234", "/path/to/EQSANS_1235.nxs.h5")
    assert configuration["DIRECT_BEAM_RUNS"] == ("EQSANS_1236",)


def test_configuration_uses_defaults_for_omitted_settings(tmp_path):
    configuration_file = tmp_path / "sensitivity.toml"
    configuration_file.write_text(
        """
FLOOD_RUNS = 1234
DIRECT_BEAM_RUNS = [1235]
SCALE_COMPONENTS = { detector1 = [1.0, 1.1, 1.0] }
"""
    )

    configuration = _SCRIPT.load_configuration(configuration_file)

    assert configuration["FLOOD_RUNS"] == (1234,)
    assert configuration["DIRECT_BEAM_RUNS"] == (1235,)
    assert configuration["DARK_CURRENT_RUNS"] is None
    assert configuration["SOLID_ANGLE_CORRECTION"] is True
    assert configuration["MIN_THRESHOLD"] == 0.5
    assert configuration["SCALE_COMPONENTS"] == {"detector1": [1.0, 1.1, 1.0]}


def test_configuration_validates_auxiliary_run_lengths():
    with pytest.raises(ValueError, match="BLOCKED_BEAM_RUNS"):
        _SCRIPT.merge_configuration(
            {
                "FLOOD_RUNS": [1, 2],
                "BLOCKED_BEAM_RUNS": [3, 4, 5],
            }
        )


def test_configuration_requires_transmission_settings_together():
    with pytest.raises(ValueError, match="must be supplied together"):
        _SCRIPT.merge_configuration({"TRANSMISSION_FLOOD_RUNS": 1})


def test_configuration_validates_scale_components():
    with pytest.raises(ValueError, match="three-number lists"):
        _SCRIPT.merge_configuration({"SCALE_COMPONENTS": {"detector1": [1.0, 1.0]}})

    with pytest.raises(ValueError, match="positive"):
        _SCRIPT.merge_configuration({"SCALE_COMPONENTS": {"detector1": [1.0, 0.0, 1.0]}})


def test_configuration_rejects_unknown_options_including_output_path():
    with pytest.raises(ValueError, match="Unknown configuration option"):
        _SCRIPT.merge_configuration({"OUTPUT_SENSITIVITY_FILE": "/tmp/sensitivity.nxs"})


def test_configuration_rejects_invalid_option_types():
    with pytest.raises(ValueError, match="SOLID_ANGLE_CORRECTION must be a boolean"):
        _SCRIPT.merge_configuration({"SOLID_ANGLE_CORRECTION": 1})

    with pytest.raises(ValueError, match="MIN_THRESHOLD cannot exceed"):
        _SCRIPT.merge_configuration({"MIN_THRESHOLD": 2.0, "MAX_THRESHOLD": 1.0})


def test_help_includes_documentation_url(capsys):
    with pytest.raises(SystemExit) as exception:
        _SCRIPT.parse_arguments(["--help"])

    assert exception.value.code == 0
    assert _SCRIPT.DOCUMENTATION_URL in capsys.readouterr().out


def test_command_line_requires_configuration_and_output():
    with pytest.raises(SystemExit):
        _SCRIPT.parse_arguments([])

    arguments = _SCRIPT.parse_arguments(["input.toml", "/tmp/sensitivity.nxs"])

    assert arguments.configuration == Path("input.toml")
    assert arguments.output == Path("/tmp/sensitivity.nxs")


class _FakePreparer:
    def __init__(self, component):
        self.component = component
        self.calls = []
        self.scale_components = None

    def set_flood_runs(self, runs):
        self.calls.append(("set_flood_runs", runs))

    def set_direct_beam_runs(self, runs):
        self.calls.append(("set_direct_beam_runs", runs))

    def set_masks(self, universal_mask, masked_pixels):
        self.calls.append(("set_masks", universal_mask, masked_pixels))

    def set_beam_center_radius(self, radius):
        self.calls.append(("set_beam_center_radius", radius))

    def set_pixel_calibration_flag(self, flag):
        self.calls.append(("set_pixel_calibration_flag", flag))

    def set_solid_angle_correction_flag(self, flag):
        self.calls.append(("set_solid_angle_correction_flag", flag))

    def set_sample_detector_offsets(self, sample_offset, detector_offset):
        self.calls.append(("set_sample_detector_offsets", sample_offset, detector_offset))

    def set_tof_cuts(self, low_tof_clip, high_tof_clip):
        self.calls.append(("set_tof_cuts", low_tof_clip, high_tof_clip))

    def set_dark_current_runs(self, runs):
        self.calls.append(("set_dark_current_runs", runs))

    def set_blocked_beam_runs(self, runs):
        self.calls.append(("set_blocked_beam_runs", runs))

    def set_empty_cell_runs(self, runs):
        self.calls.append(("set_empty_cell_runs", runs))

    def set_transmission_correction(self, flood_runs, empty_cell_runs, reference_runs):
        self.calls.append(
            (
                "set_transmission_correction",
                flood_runs,
                empty_cell_runs,
                reference_runs,
            )
        )

    def execute(self, moving_detectors, min_threshold, max_threshold, output):
        self.calls.append(("execute", moving_detectors, min_threshold, max_threshold, output))


def test_main_maps_configuration_to_preparer(monkeypatch, tmp_path):
    configuration_file = tmp_path / "sensitivity.toml"
    configuration_file.write_text(
        """
FLOOD_RUNS = [1234, "EQSANS_1235"]
DIRECT_BEAM_RUNS = "1236"
DARK_CURRENT_RUNS = 1237
BLOCKED_BEAM_RUNS = [1238, 1239]
EMPTY_CELL_RUNS = 1240
TRANSMISSION_FLOOD_RUNS = 1241
TRANSMISSION_EMPTY_CELL_RUNS = 1242
TRANSMISSION_REFERENCE_RUNS = 1243
SAMPLE_OFFSET = 285.0
DETECTOR_OFFSET = 66.714
LOW_TOF_CLIP = 1650.0
HIGH_TOF_CLIP = 3150.0
UNIVERSAL_MASK = "/path/to/mask.xml"
MASKED_PIXELS = "1-15,242-256"
MASK_BEAM_CENTER_RADIUS = 40.0
PIXEL_CALIBRATION = true
SCALE_COMPONENTS = { detector1 = [1.0, 1.1, 1.0] }
SOLID_ANGLE_CORRECTION = false
MOVING_DETECTORS = true
MIN_THRESHOLD = 0.1
MAX_THRESHOLD = 3.0
"""
    )
    preparer = _FakePreparer(component=None)

    def create_preparer(*, component):
        preparer.component = component
        return preparer

    monkeypatch.setattr(_SCRIPT, "PrepareSensitivityCorrection", create_preparer)

    _SCRIPT.main([str(configuration_file), str(tmp_path / "sensitivity.nxs")])

    assert preparer.component == "detector1"
    assert preparer.scale_components == {"detector1": [1.0, 1.1, 1.0]}
    assert preparer.calls == [
        ("set_flood_runs", (1234, "EQSANS_1235")),
        ("set_direct_beam_runs", (1236,)),
        ("set_masks", "/path/to/mask.xml", "1-15,242-256"),
        ("set_beam_center_radius", 40.0),
        ("set_pixel_calibration_flag", True),
        ("set_solid_angle_correction_flag", False),
        ("set_sample_detector_offsets", 285.0, 66.714),
        ("set_tof_cuts", 1650.0, 3150.0),
        ("set_dark_current_runs", (1237,)),
        ("set_blocked_beam_runs", (1238, 1239)),
        ("set_empty_cell_runs", (1240,)),
        ("set_transmission_correction", (1241,), (1242,), (1243,)),
        ("execute", True, 0.1, 3.0, str(tmp_path / "sensitivity.nxs")),
    ]


def test_main_uses_defaults_and_leaves_optional_corrections_disabled(monkeypatch, tmp_path):
    configuration_file = tmp_path / "sensitivity.toml"
    configuration_file.write_text(
        """
FLOOD_RUNS = 1234
DIRECT_BEAM_RUNS = 1235
"""
    )
    preparer = _FakePreparer(component=None)

    def create_preparer(*, component):
        preparer.component = component
        return preparer

    monkeypatch.setattr(_SCRIPT, "PrepareSensitivityCorrection", create_preparer)

    output_file = tmp_path / "sensitivity.nxs"
    _SCRIPT.main([str(configuration_file), str(output_file)])

    assert preparer.calls == [
        ("set_flood_runs", (1234,)),
        ("set_direct_beam_runs", (1235,)),
        ("set_masks", None, None),
        ("set_beam_center_radius", 65.0),
        ("set_pixel_calibration_flag", False),
        ("set_solid_angle_correction_flag", True),
        ("set_sample_detector_offsets", 0.0, 0.0),
        ("set_tof_cuts", 500.0, 2000.0),
        ("execute", False, 0.5, 2.0, str(output_file)),
    ]
