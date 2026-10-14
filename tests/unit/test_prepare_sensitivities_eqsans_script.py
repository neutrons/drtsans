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
