#!/usr/bin/env python
import importlib.resources

import pytest

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/mono/gpsans/attenuation.py
from drtsans.mono.gpsans import attenuation_factor
from drtsans.mono.gpsans.attenuation import (
    _ATTENUATOR_NAMES,
    _NO_ATTENUATION,
    _attenuation_factor,
    _attenuator_name,
    _load_attenuation_coefficients,
)

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py
from drtsans.samplelogs import SampleLogs


def test_attenuation_factor(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # Test input and expected values provided by Lisa Debeer-Schmitt, 2020-02-26
    wavelength = 4.75
    attenuator = 6  # x2k
    expected_value = 0.001037270673420313
    expected_error = 7.005200329552345e-05

    # Add sample logs
    SampleLogs(ws).insert("wavelength", wavelength, "Angstrom")
    SampleLogs(ws).insert("attenuator", attenuator)

    value, error = attenuation_factor(ws)
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def test_attenuation_factor_open_close(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # add wavelength
    SampleLogs(ws).insert("wavelength", 1.54, "Angstrom")

    # add Undefined attenuator
    attenuator = 0  # Undefined
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)

    # add Close attenuator
    attenuator = 1  # Close
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)

    # add Open attenuator
    attenuator = 2  # Open
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)


def test_attenuation_factor_missing_logs(generic_workspace, clean_workspace):
    """This test that correct error messages are return if the required
    attenuator or wavelenght logs is missing
    """
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # Missing attenuator and wavelength log
    with pytest.raises(RuntimeError) as excinfo:
        attenuation_factor(ws)
    assert "attenuator" in str(excinfo.value)  # Should complain about missing attenuator

    # Add in attenuator so only missing wavelength log
    SampleLogs(ws).insert("attenuator", 4)
    with pytest.raises(RuntimeError) as excinfo:
        attenuation_factor(ws)
    assert "wavelength" in str(excinfo.value)  # Should complain about missing wavelength


# Coefficients (A, A error, B, B error, C, C error) provided by Lisa Debeer-Schmitt, 2020-02-26
REFERENCE_COEFFICIENTS = {
    "x3": (
        0.3733459538730628,
        0.008609717163831113,
        0.08056544906925872,
        0.008241433507695071,
        0.0724341919138054,
        0.01125160779959418,
    ),
    "x30": (
        0.11696573514650677,
        0.006304060228295941,
        0.25014801934427583,
        0.01012469612884642,
        0.003696051816711061,
        0.0003197928933191539,
    ),
    "x300": (
        0.028719247985112162,
        0.0019190523738874328,
        0.3884993528348815,
        0.010600714703684273,
        0.00017081815634872129,
        9.055642884664314e-06,
    ),
    "x2k": (
        0.015510737042113254,
        0.0008301527045697745,
        0.5840982399579384,
        0.010252064767405953,
        6.966839283167031e-05,
        2.260503164358648e-06,
    ),
    "x10k": (
        0.00563013075327734,
        0.0005203265715819975,
        0.6961581698084675,
        0.01938010154115584,
        1.3123049075167468e-05,
        1.5266828654446554e-06,
    ),
    "x100k": (
        0.1439135754790426,
        0.005573924841205431,
        0.30824770207752383,
        0.011967728090637404,
        0.006739099792400909,
        0.0007076026868930973,
    ),
}


def test_default_coefficients_file():
    default_file = importlib.resources.files("drtsans.configuration") / "GPSANS_attenuation_coefficients.txt"
    with importlib.resources.as_file(default_file) as filename:
        coefficients = _load_attenuation_coefficients(filename)
    assert coefficients == REFERENCE_COEFFICIENTS
    assert set(coefficients) == set(_ATTENUATOR_NAMES.values()) - _NO_ATTENUATION
    # the default file is used when no file is given
    assert _load_attenuation_coefficients() == REFERENCE_COEFFICIENTS


def test_attenuation_factor_custom_file(generic_workspace, clean_workspace, tmp_path):
    ws = generic_workspace
    clean_workspace(ws)
    wavelength = 4.75
    SampleLogs(ws).insert("wavelength", wavelength, "Angstrom")
    SampleLogs(ws).insert("attenuator", 6)  # x2k

    # Custom file with x2k coefficients differing from the default ones
    custom_coefficients = (0.02, 0.001, 0.5, 0.01, 1.0e-4, 2.0e-6)
    coefficients_file = tmp_path / "custom_coefficients.txt"
    coefficients_file.write_text(
        "# attenuator name, A, A error, B, B error, C, C error\n\nx2k, " + ", ".join(map(str, custom_coefficients))
    )

    value, error = attenuation_factor(ws, coefficients_file)
    expected_value, expected_error = _attenuation_factor(*custom_coefficients, wavelength)
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)
    # the result differs from the one obtained with the default coefficients
    assert value != pytest.approx(attenuation_factor(ws)[0])

    # attenuator missing from the custom file
    SampleLogs(ws).insert("attenuator", 7)  # x10k
    with pytest.raises(ValueError, match="x10k not found"):
        attenuation_factor(ws, coefficients_file)


@pytest.mark.parametrize(
    "line, message",
    [
        ("x2k,0.1,0.01,0.5,0.01,0.001", "expected 7 comma-separated fields, found 6"),
        ("x2k,0.1,0.01,0.5,0.01,0.001,0.0001,0.1", "expected 7 comma-separated fields, found 8"),
        ("x2k,0.1,0.01,half,0.01,0.001,0.0001", "could not convert string to float"),
        (",0.1,0.01,0.5,0.01,0.001,0.0001", "empty attenuator name"),
        ("x2k,nan,0.01,0.5,0.01,0.001,0.0001", "coefficients must be finite numbers"),
        ("x2k,0.1,0.01,inf,0.01,0.001,0.0001", "coefficients must be finite numbers"),
    ],
)
def test_malformed_coefficients_file(tmp_path, line, message):
    coefficients_file = tmp_path / "malformed_coefficients.txt"
    coefficients_file.write_text("# header\n" + line + "\n")
    with pytest.raises(ValueError, match="Invalid line 2") as excinfo:
        _load_attenuation_coefficients(coefficients_file)
    assert message in str(excinfo.value)


def test_repeated_attenuator_in_coefficients_file(tmp_path):
    coefficients_file = tmp_path / "repeated_coefficients.txt"
    coefficients_file.write_text("x2k,0.1,0.01,0.5,0.01,0.001,0.0001\nx2k,0.2,0.01,0.5,0.01,0.001,0.0001\n")
    with pytest.raises(ValueError, match="Invalid line 2") as excinfo:
        _load_attenuation_coefficients(coefficients_file)
    assert "attenuator x2k is repeated" in str(excinfo.value)


def test_attenuation_factor_nonexistent_file(generic_workspace, clean_workspace, tmp_path):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", 6)  # x2k
    with pytest.raises(FileNotFoundError):
        attenuation_factor(ws, tmp_path / "nonexistent.txt")
    # the file is not read when the beam is not attenuated
    SampleLogs(ws).insert("attenuator", 2)  # Open
    assert attenuation_factor(ws, tmp_path / "nonexistent.txt") == (1, 0)


def test_attenuation_factor_unknown_attenuator(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", 9)
    with pytest.raises(ValueError, match="Unknown attenuator log value 9"):
        attenuation_factor(ws)


@pytest.mark.parametrize(
    "log_value, name",
    [
        (0, "Undefined"),
        (1, "Close"),
        (2, "Open"),
        (3, "x3"),
        (4, "x30"),
        (5, "x300"),
        (6, "x2k"),
        (7, "x10k"),
        (8, "x100k"),
    ],
)
def test_attenuator_name(generic_workspace, clean_workspace, log_value, name):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("attenuator", log_value)
    assert _attenuator_name(ws) == name


@pytest.mark.parametrize("log_value", [9, 5.5, 52.9971])
def test_attenuator_name_invalid(generic_workspace, clean_workspace, log_value):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("attenuator", log_value)
    with pytest.raises(ValueError, match="Unknown attenuator log value"):
        _attenuator_name(ws)


def test_attenuator_name_negative(generic_workspace, clean_workspace):
    """A negative log value, the x2k stage position (mm) of a run converted from a SPICE file, is Undefined"""
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", -195.999101)
    assert _attenuator_name(ws) == "Undefined"
    assert attenuation_factor(ws) == (1, 0)


if __name__ == "__main__":
    pytest.main([__file__])
